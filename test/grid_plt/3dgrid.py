# import matplotlib.pyplot as plt
import numpy as np
import sys
import vispy.scene
from vispy.scene import visuals

def process_grid(file_path, grid_shape):

    with open(file_path, 'rb') as f:
        return np.frombuffer(f.read(), dtype=np.float64).reshape(grid_shape)

def create_scatter(grid, view):
    
    grid_dims = grid.shape

    r_max, r_min = np.max(grid[0]), np.min(grid[0])
    g_max, g_min = np.max(grid[1]), np.min(grid[1])
    b_max, b_min = np.max(grid[2]), np.min(grid[2])

    total_points = grid_dims[1] * grid_dims[2] * grid_dims[3]
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
                # this makes so that they become less visible
                if r_power > 0.90:
                   continue 
                else:
                    alpha_channel = (g+b)/2

                out_coords.append([x,y,z])
                out_rgba.append((r,g,b, alpha_channel))
                

    scatter = visuals.Markers()
    scatter.set_data(np.array(out_coords), face_color=out_rgba, edge_color=(0,0,0,0.1), size=10)
    # scatter.update_gl_state(depth_test=False)
    view.add(scatter)
    view.camera = vispy.scene.TurntableCamera()
    visuals.XYZAxis(parent=view.scene)

def main():
    grid_shape = (3,60,60,60)
    canvas = vispy.scene.SceneCanvas(keys='interactive', show=True, bgcolor='w')
    view = canvas.central_widget.add_view()
    grid = process_grid('./184l/grid_30_0.5_SF0/McGrid_rec.grid', grid_shape)
    grid2 = process_grid('./184l/grid_30_0.5_SF0/McGrid_lig.grid', grid_shape)
    create_scatter(grid2,view)
    create_scatter(grid, view)
    if sys.flags.interactive != 1:
            vispy.app.run()



if __name__ == "__main__":

    main()






