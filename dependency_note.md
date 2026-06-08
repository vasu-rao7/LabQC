## 05/26 

Currently, GaussianProcesses.jl has a known [issue](https://github.com/STOR-i/GaussianProcesses.jl/issues/239) as a result of an update pushed to PDMats upon upgrade from 0.11.36 -> 0.11.37 

GP correction code will not run if GaussianProcesses.jl installs this version of PDMats when the package is added. I have avoided the issue here by first installing PDMats as a weak dependency and restricting it to 0.11.36. Then, I installed GaussianProcesses at its current version. GaussianProcesses.jl offers compatibility to older versions of PDMats. 

If we run into compatibility issues with this down the road, we'll have to see if the GaussianProcesses issue has been solved.  
