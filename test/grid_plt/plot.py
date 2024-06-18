import pyvista as pv
import numpy as np


def process_grid(file_path, grid_shape):

    with open(file_path, 'rb') as f:
        arr = np.frombuffer(f.read(), dtype=np.float64).reshape(grid_shape)

    return arr

def create_scatter(grid, p):
    
    grid_dims = grid.shape

    r_max, r_min = np.max(grid[0]), np.min(grid[0])
    g_max, g_min = np.max(grid[1]), np.min(grid[1])
    b_max, b_min = np.max(grid[2]), np.min(grid[2])

    out_rgba = []
    out_coords = []

    for x in range(grid_dims[1]): 
        for y in range(grid_dims[2]):
            for z in range(grid_dims[3]):

                r = (grid[0][x][y][z] - r_min) / (r_max - r_min)
                g = (grid[1][x][y][z] - g_min) / (g_max - g_min)
                b = (grid[2][x][y][z] - b_min) / (b_max - b_min)

                # how much red dominates over the blue and green
                r_power = abs(r - ((g+b)/2)) / r if r > 0 else 0
                # the outskirts of the grid tend to have extreme r (elec) values,
                # this makes so that they are ignored. 
                if r_power > 0.90:
                   continue 
                else:
                    alpha = (g+b)/2

                out_coords.append([x,y,z])
                out_rgba.append([r,g,b,alpha])

    rgba = np.array(out_rgba)
    coords = np.array(out_coords, dtype=np.float32)

    points = pv.wrap(coords)
    p.add_points(points,
                scalars=rgba,
                rgba=True,
                point_size=15,
                render_points_as_spheres=True)


def main():

    grid_shape = (3,60,60,60)
    grid = process_grid('./5aa9/grid_30_0.5_SF0/McGrid_rec.grid', grid_shape)

    p = pv.Plotter()
    # can also pass a ndarray slice 

    # Code for drawing a box
    #
    # size = 6 
    # coords = [[x,y,z] for x in range(size) for y in range(size) for z in range(size)]
    # cube_vertex = [[a,b,c] for a in [0, size-1] for b in [0, size-1] for c in [0, size-1]]
    # cube_edges = []
    # for vtx in cube_vertex:
    #     for i in range(3):
    #         new = vtx.copy()
    #         new[i] = 0 if vtx[i] == size-1 else size-1
    #         cube_edges.append(vtx)
    #         cube_edges.append(new)
    #
    # p.add_lines(np.array(cube_edges), color='b')
    # p.add_points(np.array(coords),
    #              color="red",
    #              render_points_as_spheres=True,
    #              point_size=15,
    #              show_edges=True)


    # Code for drawing a grid
    #
    create_scatter(grid[:,::,::,::], p)
    # altering this value affects the transparency
    # but can also impact performance
    p.enable_depth_peeling(60)

    # Code for saving a gif
    #
    # p.open_gif('lig.gif')
    # steps = 50 
    # rot = 360/steps
    # p.camera.elevation = -20
    # for _ in range(steps):
    #     p.camera.azimuth += rot 
    #     p.write_frame()
    #
    # p.close()
   
    p.show()

if __name__ == "__main__":

    main()

