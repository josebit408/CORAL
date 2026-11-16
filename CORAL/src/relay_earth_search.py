"""Search the supplied NASADEM tile for fixed-length lunar terrain-profile analogs.
Requires GDAL and NumPy. Rankings are exploratory geometry matches, not field validation.
"""
from pathlib import Path
import argparse,json,math,hashlib
import numpy as np
from relay_moon_comparison import load_gdal, require

def sampled_relay_network(distance, elevation, radius):
    """Fewest hops among sampled positions in increasing profile order.
    All interior antenna candidates are 5m high; endpoints are 2m high.
    Tie-break on largest minimum sampled interior clearance. This is not
    a 2D network optimum or an uncertainty-aware deployment design.
    """
    from relay_moon_comparison import validate_profile
    d, z = validate_profile(distance, elevation)
    require(math.isfinite(radius) and radius > 0, 'Planetary radius must be positive and finite.')
    h=np.full(len(z),5.);h[[0,-1]]=2.
    edges={}
    for i in range(len(z)-1):
        for j in range(i+1,len(z)):
            q=d[i:j+1]-d[i];length=d[j]-d[i]
            clearance=(z[i]+h[i])+(z[j]+h[j]-z[i]-h[i])*q/length-z[i:j+1]-q*(length-q)/(2*radius)
            if clearance.min()>0:edges[i,j]=float(clearance[1:-1].min()) if len(clearance)>2 else float('inf')
    hops=np.full(len(z),len(z)+1);hops[0]=0;quality=np.full(len(z),-np.inf);quality[0]=np.inf;pred={}
    for j in range(1,len(z)):
        for i in range(j):
            if (i,j) in edges:
                hh=hops[i]+1;qq=min(quality[i],edges[i,j])
                if hh<hops[j] or (hh==hops[j] and qq>quality[j]):hops[j]=hh;quality[j]=qq;pred[j]=i
    indices=[len(z)-1]
    while indices[-1]!=0:indices.append(pred[indices[-1]])
    indices=indices[::-1]
    # Adjacent positions have no interior samples. Infinity is only a neutral
    # tie-break sentinel, never a reported interior terrain clearance.
    links = [{'start_index':i, 'end_index':j,
              'min_sampled_interior_clearance_m':edges[i,j] if j > i+1 else None}
             for i,j in zip(indices,indices[1:])]
    checked = [edge['min_sampled_interior_clearance_m'] for edge in links
               if edge['min_sampled_interior_clearance_m'] is not None]
    return {'sample_indices':indices,'relay_count':len(indices)-2,
            'distance_m':[float(d[i]) for i in indices],
            'min_sampled_interior_clearance_m':min(checked) if checked else None,
            'links_without_interior_samples':sum(edge['end_index'] == edge['start_index']+1 for edge in links),
            'links':links}


def run(input_dir, moon_results, output_dir):
    """Run the fixed-tile search after checking required input files."""
    inp, moon_results, out = Path(input_dir), Path(moon_results), Path(output_dir)
    required = [inp / name for name in
                ('n31w107.hgt', 'n31w107.num', 'n31w107.swb', 'Kilbourne_profile01.gpkg')]
    required += [moon_results / name for name in ('CORAL_comparison.json', 'Site04_30m_mean.tif')]
    missing = [str(path) for path in required if not path.is_file()]
    require(not missing, 'Missing input files: ' + ', '.join(missing)
            + '. See DATA_SOURCES.md and SETUP.md for acquisition and setup.')
    old=json.loads((moon_results/'CORAL_comparison.json').read_text(encoding='utf-8-sig'))
    L=float(old['earth_metrics']['length_m'])
    require(math.isfinite(L) and L > 0, 'Reference corridor length must be positive and finite.')
    d=np.linspace(0,L,int(np.ceil(L/30))+1);t=d/L
    gdal, ogr, osr = load_gdal()
    md=gdal.Open(str(moon_results/'Site04_30m_mean.tif')); ma=md.ReadAsArray();mgt=md.GetGeoTransform()
    require(md.RasterCount == 1 and ma.ndim == 2 and np.isfinite(ma).all(),
            'The lunar mean surface must contain one finite raster band.')
    require(np.allclose(mgt, (-9000, 30, 0, 1000, 0, -30), rtol=0, atol=1e-8),
            'Expected the bundled Site04 30 m lunar mean grid and geotransform.')
    lunar_srs = osr.SpatialReference(wkt=md.GetProjection())
    require(lunar_srs.IsProjected() and abs(lunar_srs.GetSemiMajor() - 1737400) < 0.01
            and abs(lunar_srs.GetLinearUnits() - 1) < 1e-12,
            'The lunar mean surface must retain its projected Moon CRS in meters.')
    def moonline(y0):
        u=(-4000+d-mgt[0])/30-.5;v=np.full_like(d,(mgt[3]-y0)/30-.5);ii=np.floor(u).astype(int);jj=np.floor(v).astype(int);fx=u-ii;fy=v-jj
        require(((ii >= 0) & (jj >= 0) & (ii < ma.shape[1]-1) & (jj < ma.shape[0]-1)).all(),
                'Fixed lunar corridor extends beyond the supplied mean surface.')
        return (1-fx)*(1-fy)*ma[jj,ii]+fx*(1-fy)*ma[jj,ii+1]+(1-fx)*fy*ma[jj+1,ii]+fx*fy*ma[jj+1,ii+1]
    mz=moonline(-6500.);target=mz-mz[0]
    raw=gdal.Open(str(inp/'n31w107.hgt'));arr=raw.ReadAsArray()
    require(arr.ndim == 2 and min(arr.shape) >= 3, 'Earth raster must be a two-dimensional elevation grid.')
    for name in ('n31w107.num', 'n31w107.swb'):
        require((inp/name).stat().st_size == arr.size,
                f'{name}: byte count does not match the Earth elevation grid.')
    n=np.fromfile(inp/'n31w107.num',dtype=np.uint8).reshape(arr.shape);w=np.fromfile(inp/'n31w107.swb',dtype=np.uint8).reshape(arr.shape)
    good=(arr!=-32768)&(w==0)&(n>=1)&(n<=23)
    # A one-native-cell erosion rejects neighborhoods next to water, fill,
    # alternative-source pixels and edges. NUM is not a vertical-error raster.
    eroded=np.zeros_like(good)
    eroded[1:-1,1:-1]=np.logical_and.reduce([good[dy:dy+good.shape[0]-2,dx:dx+good.shape[1]-2] for dy in range(3) for dx in range(3)])
    mem=gdal.GetDriverByName('MEM').Create('',raw.RasterXSize,raw.RasterYSize,1,gdal.GDT_Byte)
    mem.SetGeoTransform(raw.GetGeoTransform());mem.SetProjection(raw.GetProjection());mem.GetRasterBand(1).WriteArray(eroded.astype('u1'))
    require(good.any(), 'No land cells pass the NASADEM source-quality screen.')
    out.mkdir(parents=True,exist_ok=True)
    demfile=out/'Earth_tile_UTM30m.tif'
    warp=gdal.Warp(str(demfile),raw,dstSRS='EPSG:32613',xRes=30,yRes=30,targetAlignedPixels=True,resampleAlg='bilinear',srcNodata=-32768,dstNodata=-9999,outputType=gdal.GDT_Float32,creationOptions=['COMPRESS=DEFLATE','TILED=YES'])
    gt=warp.GetGeoTransform();a=warp.ReadAsArray();ny,nx=a.shape
    bounds=(gt[0],gt[3]-ny*30,gt[0]+nx*30,gt[3])
    qa=gdal.Warp('',mem,format='MEM',dstSRS='EPSG:32613',outputBounds=bounds,width=nx,height=ny,resampleAlg='near',dstNodata=0).ReadAsArray().astype(bool)
    qa &= np.isfinite(a)&(a!=-9999)
    warp=None
    print(f'Prepared {nx} x {ny} UTM raster; {len(d)} equal-distance samples per {L:.3f} m profile.',flush=True)
    def sample(xs,ys):
        u=(xs-gt[0])/30-.5;v=(gt[3]-ys)/30-.5;i=np.floor(u).astype(np.int32);j=np.floor(v).astype(np.int32)
        inside=(i>=0)&(j>=0)&(i<nx-1)&(j<ny-1);ii=np.clip(i,0,nx-2);jj=np.clip(j,0,ny-2)
        valid=inside & qa[jj,ii] & qa[jj,ii+1] & qa[jj+1,ii] & qa[jj+1,ii+1]
        fx=u-i;fy=v-j
        z=(1-fx)*(1-fy)*a[jj,ii]+fx*(1-fy)*a[jj,ii+1]+(1-fx)*fy*a[jj+1,ii]+fx*fy*a[jj+1,ii+1]
        return z,valid
    def evaluate(starts,bearings):
        rad=np.radians(bearings);xs=starts[:,0,None]+np.sin(rad)[:,None]*d;ys=starts[:,1,None]+np.cos(rad)[:,None]*d
        z,valid=sample(xs,ys);valid=valid.all(axis=1)
        score=np.sqrt(np.mean(((z-z[:,0,None])-target)**2,axis=1));score[~valid]=np.inf
        return score
    def record(starts,bearings,scores):
        return np.column_stack([starts,bearings,scores])
    xx=np.arange(math.ceil((gt[0]+60)/500)*500,gt[0]+nx*30-60,500.)
    yy=np.arange(math.ceil((gt[3]-ny*30+60)/500)*500,gt[3]-60,500.)
    gx,gy=np.meshgrid(xx,yy);seeds=np.column_stack([gx.ravel(),gy.ravel()])
    records=[];coarse_total=0;coarse_valid=0
    for angle in np.arange(0,360,15.):
        for start in range(0,len(seeds),4096):
            st=seeds[start:start+4096];angles=np.full(len(st),angle);scores=evaluate(st,angles)
            coarse_total+=len(st);coarse_valid+=int(np.isfinite(scores).sum())
            ix=np.flatnonzero(np.isfinite(scores));keep=ix[np.argsort(scores[ix])[:150]]
            records.append(record(st[keep],angles[keep],scores[keep]))
        print(f'Coarse bearing {angle:.0f} degrees complete.',flush=True)
    allcoarse=np.vstack(records);allcoarse=allcoarse[np.argsort(allcoarse[:,3])]
    chosen=[]
    for r in allcoarse:
        mid=r[:2]+L/2*np.array([math.sin(math.radians(r[2])),math.cos(math.radians(r[2]))])
        if all(np.linalg.norm(mid-oldmid)>=1000 for _,oldmid in chosen):chosen.append((r,mid))
        if len(chosen)==40:break
    offsets=np.arange(-240,241,30.);ox,oy,oa=np.meshgrid(offsets,offsets,np.arange(-7.5,7.5001,1.5),indexing='ij')
    delta=np.column_stack([ox.ravel(),oy.ravel()]);ba=oa.ravel()
    fine=[];fine_total=0;fine_valid=0
    for r,mid in chosen:
        starts=r[:2]+delta;bearing=(r[2]+ba)%360;scores=evaluate(starts,bearing)
        fine_total+=len(starts);fine_valid+=int(np.isfinite(scores).sum());ix=np.flatnonzero(np.isfinite(scores));keep=ix[np.argsort(scores[ix])[:100]]
        fine.append(record(starts[keep],bearing[keep],scores[keep]))
    require(bool(fine), 'No valid coarse candidates passed the quality screen.')
    fine=np.vstack(fine);fine=fine[np.argsort(fine[:,3])]
    # Retain three geographically separated patches, not three near-duplicate
    # traces of the same ridge. Distances are between profile midpoints.
    final=[]
    for r in fine:
        mid=r[:2]+L/2*np.array([math.sin(math.radians(r[2])),math.cos(math.radians(r[2]))])
        if all(np.linalg.norm(mid-oldmid)>=5000 for _,oldmid in final):final.append((r,mid))
        if len(final)==3:break
    require(len(final)==3, 'Fewer than three separated candidates pass the search quality screen.')
    utm=osr.SpatialReference();utm.ImportFromEPSG(32613);utm.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    geo=osr.SpatialReference();geo.ImportFromEPSG(4326);geo.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    togeo=osr.CoordinateTransformation(utm,geo)
    def metrics(z):
        direct=(z[0]+2)*(1-t)+(z[-1]+2)*t;bulge=d*(L-d)/(2*6371000.);obs=z+bulge-direct
        j=int(np.argmax(z));m={'relief_m':float(np.ptp(z)),'max_obstruction_m':float(obs.max()),'max_obstruction_distance_m':float(d[np.argmax(obs)]),'peak_index':j,'peak_distance_m':float(d[j]),'rmse_m':float(np.sqrt(np.mean(((z-z[0])-target)**2))),'neighbor_reference_rmse_m':{str(int(y)):float(np.sqrt(np.mean(((z-z[0])-(moonline(y)-moonline(y)[0]))**2))) for y in [-6000.,-6500.,-7000.]}}
        if 0<j<len(d)-1:
            p=d[j];rt=z[j]+5
            l=(z[0]+2)+(rt-z[0]-2)*d[:j+1]/p-z[:j+1]-d[:j+1]*(p-d[:j+1])/(2*6371000.)
            dr=d[j:]-p;rr=rt+(z[-1]+2-rt)*dr/(L-p)-z[j:]-dr*(L-p-dr)/(2*6371000.)
            m.update({'peak_relay_height_m':5.,'peak_relay_min_clearance_left_m':float(l.min()),'peak_relay_min_clearance_right_m':float(rr.min()),'peak_relay_both_legs_clear':bool(l.min()>0 and rr.min()>0)})
        return m,direct,obs
    results=[]
    def native_qa(xs,ys):
        ll=np.array(togeo.TransformPoints(list(zip(map(float,xs),map(float,ys)))))
        rg=raw.GetGeoTransform();u=(ll[:,0]-rg[0])/rg[1]-.5;v=(ll[:,1]-rg[3])/rg[5]-.5;ii=np.floor(u).astype(int);jj=np.floor(v).astype(int)
        require((ii>=1).all() and (jj>=1).all() and (ii<arr.shape[1]-1).all() and (jj<arr.shape[0]-1).all(),
                'A selected candidate extends beyond native NASADEM support bounds.')
        vals=np.concatenate([n[jj,ii],n[jj,ii+1],n[jj+1,ii],n[jj+1,ii+1]])
        waters=np.concatenate([w[jj,ii],w[jj,ii+1],w[jj+1,ii],w[jj+1,ii+1]])
        return {'source_num_values':np.unique(vals).tolist(),'all_srtm_support':bool(((vals>=1)&(vals<=23)).all()),'water_support_pixels':int((waters!=0).sum())}
    for idx,(r,mid) in enumerate(final,1):
        x0,y0,bearing,score=map(float,r);rad=math.radians(bearing);xs=x0+d*math.sin(rad);ys=y0+d*math.cos(rad);z,valid=sample(xs,ys);require(valid.all(), 'Selected candidate failed grid-quality checks.')
        m,direct,obs=metrics(z);require(abs(m['rmse_m']-score)<1e-6, 'Candidate score consistency check failed.')
        latlon=togeo.TransformPoint(float(mid[0]),float(mid[1]));qaresult=native_qa(xs,ys);require(qaresult['all_srtm_support'] and qaresult['water_support_pixels']==0, 'Selected candidate failed native NASADEM quality checks.')
        m.update({'candidate_id':f'RS{idx:02d}','bearing_deg_clockwise_from_north':bearing,'start_xy_m':[x0,y0],'end_xy_m':[float(xs[-1]),float(ys[-1])],'midpoint_lon_lat':[latlon[0],latlon[1]],'qa':qaresult})
        # Local diagnostic around the selected line. These are still the same
        # DEM and reference data, not independent validation or error estimates.
        diagnostic=[]
        for offset in (-60.,0.,60.):
            for angle_delta in (-1.5,0.,1.5):
                st=np.array([[x0+offset*math.cos(rad),y0-offset*math.sin(rad)]])
                sc=evaluate(st,np.array([(bearing+angle_delta)%360]))[0]
                if np.isfinite(sc):diagnostic.append(float(sc))
        m['local_position_bearing_diagnostic_rmse_range_m']=[min(diagnostic),max(diagnostic)]
        m['local_position_bearing_diagnostic_valid_count']=len(diagnostic)
        results.append({'metrics':m,'x_m':xs.tolist(),'y_m':ys.tolist(),'elevation_m':z.tolist(),'direct_line_m':direct.tolist(),'obstruction_m':obs.tolist(),'nominal_sampled_relay_network':sampled_relay_network(d,z,6371000.)})
    # Re-evaluate the ORIGINAL saved Kilbourne line on the SAME full-tile grid.
    baseline=ogr.Open(str(inp/'Kilbourne_profile01.gpkg'))
    require(baseline is not None and baseline.GetLayerCount() > 0, 'Cannot open the Kilbourne baseline layer.')
    feature=baseline.GetLayer(0).GetNextFeature()
    require(feature is not None and feature.GetGeometryRef() is not None, 'The Kilbourne baseline layer is empty.')
    geom=feature.GetGeometryRef().Clone()
    require(geom.GetGeometryType() in (ogr.wkbLineString, ogr.wkbLineString25D) and geom.GetPointCount() >= 2,
            'The Kilbourne baseline must be a line with at least two vertices.')
    sx,sy,_=geom.GetPoint(0);ex,ey,_=geom.GetPoint(geom.GetPointCount()-1)
    require(abs(math.hypot(ex-sx,ey-sy)-L)<1e-6, 'Baseline length differs from the fixed lunar reference length.')
    bx=sx+(ex-sx)*t;by=sy+(ey-sy)*t;bz,bv=sample(bx,by);require(bv.all(), 'Baseline failed grid-quality checks.');bm,bd,bo=metrics(bz)
    # Same-time samples and no vertical scaling for every score, including baseline.
    gpkg=out/'CORAL_Earth_candidates.gpkg'
    if gpkg.exists():gpkg.unlink()
    g=ogr.GetDriverByName('GPKG').CreateDataSource(str(gpkg))
    def layer(name,typ,fields):
        lay=g.CreateLayer(name,utm,typ)
        for field,kind in fields:lay.CreateField(ogr.FieldDefn(field,kind))
        return lay
    lines=layer('candidate_profiles',ogr.wkbLineString25D,[('candidate_id',ogr.OFTString),('rmse_m',ogr.OFTReal),('relief_m',ogr.OFTReal),('obstruction_m',ogr.OFTReal),('bearing_deg',ogr.OFTReal),('status',ogr.OFTString)])
    points=layer('candidate_trial_points',ogr.wkbPoint25D,[('candidate_id',ogr.OFTString),('name',ogr.OFTString),('antenna_height_m',ogr.OFTReal),('ground_elev_m',ogr.OFTReal),('status',ogr.OFTString)])
    for p in results:
        m=p['metrics'];f=ogr.Feature(lines.GetLayerDefn())
        for name,value in [('candidate_id',m['candidate_id']),('rmse_m',m['rmse_m']),('relief_m',m['relief_m']),('obstruction_m',m['max_obstruction_m']),('bearing_deg',m['bearing_deg_clockwise_from_north']),('status','Exploratory geometry match; land cover and access not validated')]:f.SetField(name,value)
        geom=ogr.Geometry(ogr.wkbLineString25D)
        for x,y,z in zip(p['x_m'],p['y_m'],p['elevation_m']):geom.AddPoint(x,y,z)
        f.SetGeometry(geom);lines.CreateFeature(f)
        for name,i,h in [('Base',0,2.),('Peak_relay_trial',m['peak_index'],5.),('Rover',len(d)-1,2.)]:
            f=ogr.Feature(points.GetLayerDefn());f.SetField('candidate_id',m['candidate_id']);f.SetField('name',name);f.SetField('antenna_height_m',h);f.SetField('ground_elev_m',p['elevation_m'][i]);f.SetField('status','Trial point, not deployment recommendation')
            geom=ogr.Geometry(ogr.wkbPoint25D);geom.AddPoint(p['x_m'][i],p['y_m'][i],p['elevation_m'][i]);f.SetGeometry(geom);points.CreateFeature(f)
    bl=layer('Kilbourne_baseline',ogr.wkbLineString25D,[('name',ogr.OFTString),('rmse_m',ogr.OFTReal),('obstruction_m',ogr.OFTReal)])
    f=ogr.Feature(bl.GetLayerDefn());f.SetField('name','Original Kilbourne trial');f.SetField('rmse_m',bm['rmse_m']);f.SetField('obstruction_m',bm['max_obstruction_m']);geom=ogr.Geometry(ogr.wkbLineString25D)
    for x,y,z in zip(bx,by,bz):geom.AddPoint(float(x),float(y),float(z))
    f.SetGeometry(geom);bl.CreateFeature(f)
    focus=layer('RS01_focus_profile',ogr.wkbLineString25D,[('candidate_id',ogr.OFTString),('rmse_m',ogr.OFTReal)])
    p=results[0];f=ogr.Feature(focus.GetLayerDefn());f.SetField('candidate_id','RS01');f.SetField('rmse_m',p['metrics']['rmse_m']);geom=ogr.Geometry(ogr.wkbLineString25D)
    for x,y,z in zip(p['x_m'],p['y_m'],p['elevation_m']):geom.AddPoint(x,y,z)
    f.SetGeometry(geom);focus.CreateFeature(f)
    antennas=layer('nominal_network_antennas',ogr.wkbPoint25D,[('candidate_id',ogr.OFTString),('name',ogr.OFTString),('mast_m',ogr.OFTReal),('ground_m',ogr.OFTReal),('top_m',ogr.OFTReal),('status',ogr.OFTString)])
    links=layer('nominal_network_links',ogr.wkbLineString25D,[('candidate_id',ogr.OFTString),('link',ogr.OFTInteger),('clearance_m',ogr.OFTReal),('status',ogr.OFTString)])
    for p in results:
        cid=p['metrics']['candidate_id'];indices=p['nominal_sampled_relay_network']['sample_indices']
        for counter,i in enumerate(indices):
            mast=2. if i in (0,len(d)-1) else 5.;top=p['elevation_m'][i]+mast
            f=ogr.Feature(antennas.GetLayerDefn());f.SetField('candidate_id',cid);f.SetField('name','Base' if i==0 else 'Rover' if i==len(d)-1 else f'R{counter:02d}');f.SetField('mast_m',mast);f.SetField('ground_m',p['elevation_m'][i]);f.SetField('top_m',top);f.SetField('status','Nominal sample-only model; elevation uncertainty and RF not modeled')
            geom=ogr.Geometry(ogr.wkbPoint25D);geom.AddPoint(p['x_m'][i],p['y_m'][i],top);f.SetGeometry(geom);antennas.CreateFeature(f)
        for counter,edge in enumerate(p['nominal_sampled_relay_network']['links'],1):
            f=ogr.Feature(links.GetLayerDefn());f.SetField('candidate_id',cid);f.SetField('link',counter)
            if edge['min_sampled_interior_clearance_m'] is not None:
                f.SetField('clearance_m',edge['min_sampled_interior_clearance_m'])
            f.SetField('status','No interior terrain samples' if edge['min_sampled_interior_clearance_m'] is None else 'Nominal sample-only geometric visibility')
            geom=ogr.Geometry(ogr.wkbLineString25D)
            for i in (edge['start_index'],edge['end_index']):geom.AddPoint(p['x_m'][i],p['y_m'][i],p['elevation_m'][i]+(2. if i in (0,len(d)-1) else 5.))
            f.SetGeometry(geom);links.CreateFeature(f)
    g=None
    qa_counts={'raw_water_pixels':int((w==255).sum()),'raw_interpolation_pixels':int((n==250).sum()),'raw_error_pixels':int((n==255).sum()),'raw_srtm_pixels':int(((n>=1)&(n<=23)).sum()),'raw_total_pixels':int(n.size)}
    result={'purpose':'Explore Earth terrain analog candidates for one fixed lunar communications-obstruction profile','input_sources':[{'file':(inp/name).name,'sha256':hashlib.sha256((inp/name).read_bytes()).hexdigest()} for name in ['n31w107.hgt','n31w107.num','n31w107.swb','Kilbourne_profile01.gpkg']],'reference_result_sha256':hashlib.sha256((moon_results/'CORAL_comparison.json').read_bytes()).hexdigest(),'length_m':L,'distance_m':d.tolist(),'sample_step_m':float(d[1]),'sample_count':len(d),'moon_elevation_m':mz.tolist(),'moon_primary_start_xy_m':[-4000.,-6500.],'moon_primary_end_xy_m':[-4000.+L,-6500.],'selection_metric':'Root mean squared difference of start-relative ground heights in meters. Lower is closer. No vertical scaling, no path stretching. Primary reference only drives ranking.','search_parameters':{'coarse_start_spacing_m':500,'coarse_bearing_step_deg':15,'fine_seed_count':len(chosen),'fine_seed_midpoint_separation_m':1000,'fine_start_offsets_m':[-240,240],'fine_start_step_m':30,'fine_bearing_offsets_deg':[-7.5,7.5],'fine_bearing_step_deg':1.5,'final_midpoint_separation_m':5000,'tile_wgs84_bounds':[-107,31,-106,32]},'counts':{'coarse_evaluated':coarse_total,'coarse_quality_valid':coarse_valid,'fine_evaluated':fine_total,'fine_quality_valid':fine_valid},'quality':qa_counts,'quality_method':'Require SRTM-origin NUM 1–23, SWB land, valid heights, one-cell native quality erosion, and all four bilinear sampling supports passing the reprojected mask. Final native support codes checked again. NUM is source/count metadata, not a height-uncertainty estimate.','sampling':'Entire Earth tile bilinearly reprojected to aligned EPSG:32613 30m cells. Both planets sampled at 75 equal-distance locations (~29.76m apart). Lunar fixed 30m mean DEM as previously documented. Original baseline re-evaluated on the same Earth grid.','physics':'2m base/rover, 5m peak-relay trial. Approximate spherical curvature, Earth R=6371000m. No atmosphere, Fresnel-zone or RF propagation model.','baseline_metrics':bm,'baseline_x_m':bx.tolist(),'baseline_y_m':by.tolist(),'baseline_elevation_m':bz.tolist(),'baseline_obstruction_m':bo.tolist(),'candidates':results,'limitations':['Geometry-only retrieval within one 1-degree tile, not a global optimum or verified new analog discovery.','Ranked matches are selected on the same data used to report fit. Nearby lunar references are secondary diagnostics, not independent validation.','Land cover, geology, human structures, accessibility, ownership and permissions remain unverified.','Water/source masks do not establish a bare-earth terrain model or an operational site.','Earth elevation uncertainty is not quantified. Lower RMSE is not a probability of mission success.','30m grid spacing does not imply 30m effective resolution or equivalent accuracy across planets.','Trial lines are antenna links, not drivable rover routes.','All marker elevations are ground heights; antenna heights are attributes.']}
    result['moon_nominal_sampled_relay_network']=sampled_relay_network(d,mz,1737400.)
    result['baseline_nominal_sampled_relay_network']=sampled_relay_network(d,bz,6371000.)
    result['network_method']='Minimum hops on an increasing-index graph of the same 75 sample positions; endpoints 2m, relay candidates 5m. Edges require positive nominal sampled clearance including curvature. Tie-break maximizes the minimum sampled interior clearance. Not a deployment design, a 2D optimum, or uncertainty-aware.'
    result['moon_max_direct_obstruction_m']=float(np.max(mz+d*(L-d)/(2*1737400.)-np.interp(d,[0,L],[mz[0]+2,mz[-1]+2])))
    result['moon_relief_m']=float(np.ptp(mz))
    (out/'CORAL_Earth_search.json').write_text(json.dumps(result,indent=2,allow_nan=False),encoding='utf-8')
    check=ogr.Open(str(gpkg));expected=[3,9,1,1,sum(len(p['nominal_sampled_relay_network']['sample_indices']) for p in results),sum(len(p['nominal_sampled_relay_network']['links']) for p in results)];require([check.GetLayer(i).GetFeatureCount() for i in range(6)]==expected, 'Exported GIS feature counts differ from the calculated networks.')
    print(json.dumps({'counts':result['counts'],'baseline':bm,'candidates':[v['metrics'] for v in results]},indent=2),flush=True)

def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir', type=Path, default=Path('data/inputs'),
                        help='NASADEM tile, masks and Kilbourne GeoPackage (default: data/inputs).')
    parser.add_argument('--moon-results', type=Path, default=Path('data/reference_moon'),
                        help='Prepared lunar reference folder (default: data/reference_moon).')
    parser.add_argument('--output-dir', type=Path, default=Path('generated/earth'),
                        help='Generated rasters and search outputs (default: generated/earth).')
    args = parser.parse_args(argv)
    try:
        run(args.input_dir, args.moon_results, args.output_dir)
    except (OSError, ValueError, RuntimeError, KeyError) as error:
        parser.exit(2, f'CORAL search error: {error}\n')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
