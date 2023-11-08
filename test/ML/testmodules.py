import grid_loader 

data_loader = grid_loader.GridLoader("/data/alex/workproj/grids", (7, 20, 20, 20))
arr = data_loader.create_grids_ndarray()

print(arr)




