# Changelog

## 2026-11-16

### Packaging and documentation

- Added the published demo and repository links to project materials.
- Updated setup instructions for the existing GitHub repository and its `/docs` Pages deployment.
- Added a Conda environment specification, minimal test requirements and a GitHub Actions workflow.
- Documented the difference between the original 83-point comparison and the current 75-point candidate search.
- Added reference-data and results-field documentation with units and interpretation limits.

### Code and demo fixes

- Made Earth analysis scripts safe to import and aligned their default paths with the repository structure.
- Added clear missing-input errors and checks for lunar grids, source masks and GIS baseline inputs.
- Corrected unchecked-link metadata: a link between adjacent samples has `null` interior clearance. Added counts of links without interior samples. Updated the published JSON and GIS link attribute consistently.
- Preserved all published elevation arrays, candidate scores, relay counts and sampled relay positions.
- Improved demo keyboard focus, navigation, small-screen chart access, source attribution and JSON downloads.
- Extended tests to compare the demo's embedded profiles and metrics with the published results, and to check CLI behavior and relay edge cases.

### Validation scope

The numerical suite passes locally, including 36 demo candidate/height/view combinations checked against a separate Python calculation. GIS hillshade preparation and figure generation were exercised with existing processed data. The full Earth search was not rerun for this packaging update. Browser visual inspection and independent elevation/field validation remain outstanding. The GitHub Actions workflow will run after upload and has not yet run on the remote repository.
