"""Render search results. Requires NumPy, Matplotlib and Pillow."""
from pathlib import Path
import json,argparse
import numpy as np
from PIL import Image
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ap=argparse.ArgumentParser();ap.add_argument('results_dir',nargs='?',default='results');args=ap.parse_args();out=Path(args.results_dir)
r=json.loads((out/'RelayScout_Earth_search.json').read_text());d=np.array(r['distance_m']);moon=np.array(r['moon_elevation_m']);cands=r['candidates'];colors=['#16866d','#326db1','#b47724']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11})

# Exact same physical axes and vertical scale for all profiles.
fig,(ax,bx)=plt.subplots(2,1,figsize=(11,8.2),gridspec_kw={'height_ratios':[1.7,1]})
fig.subplots_adjust(left=.1,right=.97,top=.84,bottom=.20,hspace=.42)
fig.suptitle('RelayScout | Earth terrain candidate search',x=.055,y=.98,ha='left',fontsize=18,fontweight='bold')
fig.text(.055,.925,'Fixed lunar reference. Same 2.20 km length and 75 sample positions. No vertical stretching.',fontsize=11,color='#4c5967')
ax.plot(d,moon-moon[0],color='#202c39',linewidth=2.7,label='Moon reference',zorder=6)
for p,color in zip(cands,colors):
    z=np.array(p['elevation_m']);ax.plot(d,z-z[0],color=color,linewidth=1.8,label=p['metrics']['candidate_id'])
base=np.array(r['baseline_elevation_m']);ax.plot(d,base-base[0],color='#7f8790',linestyle='--',linewidth=1.8,label='Kilbourne baseline')
ax.set(ylabel='Height relative to start (m)',xlabel='Distance from trial base (m)',title='Retrieved profiles more closely resemble the reference terrain')
ax.legend(frameon=False,fontsize=9,ncol=3,loc='upper left');ax.grid(alpha=.2)
scores=[r['baseline_metrics']['rmse_m']]+[p['metrics']['rmse_m'] for p in cands]
bars=bx.barh(['Kilbourne baseline','RS01','RS02','RS03'],scores,color=['#929aa3']+colors,height=.65);bx.invert_yaxis();bx.set_xlim(0,250)
for rect,score in zip(bars,scores):bx.text(score+3,rect.get_y()+rect.get_height()/2,f'{score:.1f} m',va='center',fontsize=11)
bx.set(xlabel='Profile-shape RMSE (m). Lower is closer.',title='Fit metric used for ranking; this is not positional accuracy or a similarity percentage')
bx.grid(axis='x',alpha=.2);bx.set_axisbelow(True)
for a in (ax,bx):a.spines[['top','right']].set_visible(False)
fig.text(.055,.095,'Search domain: the supplied 1-degree Earth tile. Water and non-SRTM source neighborhoods excluded.\nCandidates are selected on these same data; independent imagery/terrain and field validation remain outstanding.',fontsize=10,color='#33475a')
fig.text(.055,.025,'NASA/USGS NASADEM and NASA LRO–LOLA. Approximately 30 m grids do not imply equal accuracy or effective resolution.\nTerrain geometry only: no inference of geological equivalence, site access, radio performance or newly discovered analog status.',fontsize=9,color='#58626b')
fig.savefig(out/'RelayScout_candidate_comparison.png',dpi=170,bbox_inches='tight');plt.close(fig)

# A real DEM-derived map, with no invented roads, political boundaries or imagery.
meta=json.loads((out/'display_grid.json').read_text());gt=meta['geotransform'];shade=np.array(Image.open(out/'Earth_tile_hillshade.tif'));h,w=shade.shape
extent=(gt[0]/1000,(gt[0]+30*w)/1000,(gt[3]-30*h)/1000,gt[3]/1000)
masked=np.ma.masked_equal(shade,0);cmap=plt.get_cmap('gray').copy();cmap.set_bad('white')
fig,axes=plt.subplots(2,2,figsize=(11,10));fig.subplots_adjust(left=.08,right=.97,top=.88,bottom=.12,wspace=.26,hspace=.32)
fig.suptitle('RelayScout | Candidate locations',x=.055,y=.98,ha='left',fontsize=18,fontweight='bold')
fig.text(.055,.932,'UTM zone 13N. Line lengths approximately 2.20 km. Hillshade is simulated lighting for display.',fontsize=10,color='#4c5967')
whole=axes.flat[0];whole.imshow(masked,cmap=cmap,extent=extent,origin='upper',vmin=1,vmax=255)
for p,color in zip(cands,colors):
    x=np.array(p['x_m'])/1000;y=np.array(p['y_m'])/1000;whole.plot(x,y,color=color,linewidth=2.3);whole.scatter([x.mean()],[y.mean()],s=25,color=color);whole.annotate(p['metrics']['candidate_id'],(x.mean(),y.mean()),xytext=(9,3),textcoords='offset points',fontsize=10,color=color,fontweight='bold')
whole.plot(np.array(r['baseline_x_m'])/1000,np.array(r['baseline_y_m'])/1000,color='#743c8b',linewidth=2);whole.annotate('Kilbourne',(r['baseline_x_m'][0]/1000,r['baseline_y_m'][0]/1000),xytext=(4,-13),textcoords='offset points',fontsize=9,color='#743c8b')
whole.set_title('Uploaded tile and shortlist',fontsize=12)
for ax,p,color in zip(list(axes.flat)[1:],cands,colors):
    x=np.array(p['x_m'])/1000;y=np.array(p['y_m'])/1000;cx=(x[0]+x[-1])/2;cy=(y[0]+y[-1])/2
    # Crop image pixels to the displayed area to avoid repeated full-tile rasterization.
    half=2.;xmin,xmax=cx-half,cx+half;ymin,ymax=cy-half,cy+half
    i0=max(0,int((xmin*1000-gt[0])/30));i1=min(w,int((xmax*1000-gt[0])/30)+1);j0=max(0,int((gt[3]-ymax*1000)/30));j1=min(h,int((gt[3]-ymin*1000)/30)+1)
    ex=( (gt[0]+30*i0)/1000,(gt[0]+30*i1)/1000,(gt[3]-30*j1)/1000,(gt[3]-30*j0)/1000 )
    ax.imshow(masked[j0:j1,i0:i1],cmap=cmap,extent=ex,origin='upper',vmin=1,vmax=255);ax.plot(x,y,color=color,linewidth=2.4)
    ax.scatter([x[0],x[-1]],[y[0],y[-1]],s=35,c=[color,color],edgecolors='white',linewidths=.7,zorder=4)
    for name,i,offset in [('Base',0,(5,8)),('Rover',-1,(5,-15))]:ax.annotate(name,(x[i],y[i]),xytext=offset,textcoords='offset points',fontsize=9,color='#152a3b',bbox={'facecolor':'white','edgecolor':'none','alpha':.75,'pad':1})
    lon,lat=p['metrics']['midpoint_lon_lat'];ax.set_title(f"{p['metrics']['candidate_id']} | midpoint {lat:.4f}, {lon:.4f}",fontsize=11,color=color);ax.set_xlim(xmin,xmax);ax.set_ylim(ymin,ymax)
for ax in axes.flat:ax.set_xlabel('Easting (km)');ax.set_ylabel('Northing (km)');ax.ticklabel_format(useOffset=False,style='plain');ax.set_aspect('equal')
fig.text(.055,.043,'Coordinates identify computational terrain candidates. Land cover, geology, ownership and access have not been verified.\nA line is a communications test geometry, not a proposed hiking or driving route. Grid north is up.',fontsize=10,color='#4c5967')
fig.savefig(out/'RelayScout_candidate_map.png',dpi=170,bbox_inches='tight');plt.close(fig)

# Show the nominal sampled-position network for the best geometry candidate.
p=cands[0];z=np.array(p['elevation_m']);net=p['nominal_sampled_relay_network'];ids=net['sample_indices'];tops=z+5;tops[[0,-1]]=z[[0,-1]]+2
fig,ax=plt.subplots(figsize=(11,6));fig.subplots_adjust(left=.09,right=.97,top=.83,bottom=.28)
fig.suptitle('RelayScout | RS01 nominal relay demonstration',x=.055,y=.98,ha='left',fontsize=17,fontweight='bold')
fig.text(.055,.92,f"{net['relay_count']} intermediate 5 m relays among 75 tested positions. Base and rover antennas: 2 m.",fontsize=11,color='#4c5967')
ax.fill_between(d,z- z[0],-180,color='#e1e5e9');ax.plot(d,z-z[0],color='#45566b',linewidth=1.8,label='RS01 ground')
ax.plot(d,moon-moon[0],color='#8a6b49',linestyle=':',linewidth=1.5,label='Moon reference ground')
ax.plot(d,np.array(p['direct_line_m'])-z[0],color='#c14a43',linestyle='--',linewidth=1.5,label='Obstructed direct path')
ax.plot(d[ids],tops[ids]-z[0],color=colors[0],linewidth=2,label='Nominal relay links');ax.scatter(d[ids],tops[ids]-z[0],color=colors[0],s=28,zorder=5)
for counter,i in enumerate(ids):
    name='Base' if i==0 else 'Rover' if i==len(d)-1 else f'R{counter}'
    shift=(-14 if counter%2 else 14);ha='right' if i==len(d)-1 else 'left';dx=3
    if i==0:shift=22
    if counter==1:shift=-18
    if i==len(d)-1:shift=-17;dx=-3
    if counter==len(ids)-2:shift=15;dx=-10;ha='right'
    ax.annotate(name,(d[i],tops[i]-z[0]),xytext=(dx,shift),textcoords='offset points',fontsize=9,ha=ha,color='#1d4037')
ax.set(xlabel='Distance from trial base (m)',ylabel='Height relative to base ground (m)',xlim=(-30,d[-1]+30),ylim=(-165,360));ax.grid(alpha=.2);ax.spines[['top','right']].set_visible(False)
fig.legend(*ax.get_legend_handles_labels(),loc='lower left',bbox_to_anchor=(.075,.095),frameon=False,ncol=2,fontsize=10)
fig.text(.055,.025,'A discrete, single-profile geometric model, not a deployment design. Heights are vertically exaggerated.\nEarth elevation uncertainty, subgrid obstacles, Fresnel-zone clearance and radio propagation are not modeled.',fontsize=9,color='#58626b')
fig.savefig(out/'RelayScout_RS01_network.png',dpi=170,bbox_inches='tight');plt.close(fig)
print('Rendered comparison, candidate map and nominal relay demonstration.')
