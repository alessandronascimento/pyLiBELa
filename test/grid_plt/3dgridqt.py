import numpy as np
import sys
import pyqtgraph as pg 
import pyqtgraph.opengl as gl

def process_grid(file_path, grid_shape):

    with open(file_path, 'rb') as f:
        return np.frombuffer(f.read(), dtype=np.float64).reshape(grid_shape)

def create_scatter(grid, view: gl.GLViewWidget):
    
    grid_dims = grid.shape

    r_max, r_min = np.max(grid[0]), np.min(grid[0])
    g_max, g_min = np.max(grid[1]), np.min(grid[1])
    b_max, b_min = np.max(grid[2]), np.min(grid[2])

    out_rgba = []
    out_coords = []

    scatter = gl.GLScatterPlotItem()
    view.addItem(scatter)

    for x in range(grid_dims[1]): 
        for y in range(grid_dims[2]):
            for z in range(grid_dims[3]):

                r = (grid[0][x][y][z] - r_min) / (r_max - r_min)
                g = (grid[1][x][y][z] - g_min) / (g_max - g_min)
                b = (grid[2][x][y][z] - b_min) / (b_max - b_min)

                # how much red dominate over the blue and green
                r_power = abs(r - ((g+b)/2)) / r if r > 0 else 0
                # the outskirts of the grid tend to have extreme r (elec) values,
                # this makes so that they become less visible
                if r_power > 0.90:
                   continue 
                else:
                    alpha_channel = (g+b)/2

                out_coords.append([x,y,z])
                out_rgba.append([r,g,b, alpha_channel])

    out_coords = np.array(out_coords)
    out_rgba = np.array(out_rgba)
                
    scatter.setGLOptions('translucent')
    scatter.setData(pos=out_coords, color=out_rgba)


def main():

    app = pg.mkQApp()
    view = gl.GLViewWidget()
    view.setBackgroundColor(0.2)
    view.show()

    grid_shape = (3,60,60,60)
    grid = process_grid('./10gs/grid_30_0.5_SF0/McGrid_rec.grid', grid_shape)
    create_scatter(grid[:,:,:,:], view)

    sys.exit(app.exec_())


if __name__ == "__main__":

    main()






