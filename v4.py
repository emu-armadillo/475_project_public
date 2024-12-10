import laspy
import open3d as o3d
import numpy as np

import pyvista as pv
#yes
import os as os
import random


    



#splits point cloud
def one_ninth_reduction():
    reduce_point_cloud = o3d.io.read_point_cloud("regular.ply")
    bbox = reduce_point_cloud.get_axis_aligned_bounding_box()
    print(f"bounding box{bbox}")


    x_dim = bbox.max_bound[0] - bbox.min_bound[0]
    y_dim = bbox.max_bound[1] - bbox.min_bound[1]
    z_dim = bbox.max_bound[2] - bbox.min_bound[2]

    min_bound = np.array([bbox.min_bound[0]+((x_dim*0)/3),bbox.min_bound[1],bbox.min_bound[2]])
    max_bound = np.array([bbox.min_bound[0]+((x_dim*1)/3),bbox.min_bound[1]+(y_dim/3),bbox.max_bound[2]])
    cropped_pcd = reduce_point_cloud.crop(o3d.geometry.AxisAlignedBoundingBox(min_bound, max_bound))
    o3d.visualization.draw_geometries( [cropped_pcd])
    o3d.io.write_point_cloud("regular5.ply",cropped_pcd)


def pv_decimate_meshv2(infile, outfile):
    las = laspy.read(infile)
    
    print("read done")
    # Create a PyVista PolyData object from the points
    point_cloud = pv.PolyData(las.xyz)

    #print("conversion1 done")

    tolerance = 0.02
    #cloud.plot(point_size=2)
    cloud2 = point_cloud.clean(
    point_merging=True,
    merge_tol=tolerance,
    lines_to_points=False,
    polys_to_lines=False,
    strips_to_polys=False,
    inplace=False,
    absolute=False,
    progress_bar=True,  
    )

  

    # surf = point_cloud.delaunay_2d()
    # point_cloud = None
    # print("conversion 2 done")
    # decimated = surf.decimate(0.999)
    # print("decimation done")
    # # Visualize the point cloud
    # #point_cloud.plot(style='points_gaussian', point_size=2) 
    # p = pv.Plotter()
    # #p.add_mesh(point_cloud)
    # p.add_mesh(decimated)
    # p.show()
    surf = cloud2.delaunay_2d()
    #decimated = surf.decimate(0.999)

    p = pv.Plotter()
    p.add_mesh(surf)
    p.add_mesh(point_cloud,point_size = 1, color ="red")
    p.show()
    cloud2.save(outfile)




def pv_decimate_meshv3(infile, outfile):
    cloud = pv.read("r_test_points.ply")
    print(f"number of points {cloud.number_of_points}")
#     tolerance = 0.01
#    # cloud.plot(point_size=2)
#     cloud2 = cloud.clean(
#     point_merging=True,
#     merge_tol=tolerance,
#     lines_to_points=False,
#     polys_to_lines=False,
#     strips_to_polys=False,
#     inplace=False,
#     absolute=False,
#     progress_bar=True,  
#     )

    #print(f"number of points after reduction {cloud2.number_of_points}")

    print("loaded")
    surf = cloud.delaunay_2d(progress_bar=True)

    reduced1 = surf.decimate(0.99,progress_bar=True)
    #decimated = surf.decimate(0.999)

    # p = pv.Plotter()
    # p.add_mesh(cloud,color='red',point_size=1)
    # p.add_mesh(surf)
    # p.show()
    # p = pv.Plotter()
    # p.add_mesh(cloud,color='red',point_size=1)
    # p.add_mesh(reduced1)
    # p.show()

    num_ver = len(reduced1.points)
    print(f"number of vertexs in the mesh is {num_ver}")

    save_cloud =  pv.PolyData(reduced1.points)
    pv.PolyData.save(save_cloud, "r_test_points.ply")





#pyvist decimate mesh
def pv_decimate_mesh(infile):
      #las = laspy.read(infile)
    outfile = infile.replace(".ply","_reduced.ply")
    
    print("read done")
    # Create a PyVista PolyData object from the points
    cloud = pv.read(infile)
    scale_factor = 4
    #this is to ensure the plane is created in the correct direction
    cloud = cloud.scale([scale_factor,scale_factor,1])


   
    surf = cloud.delaunay_2d(progress_bar=True)

    reduced1 = surf.decimate(0.99,progress_bar=True)
    #decimated = surf.decimate(0.999)

    num_ver = len(reduced1.points)


    cloud = cloud.scale([1/scale_factor,1/scale_factor,1])
    reduced1 = reduced1.scale([1/scale_factor,1/scale_factor,1])



    save_cloud =  pv.PolyData(reduced1.points)
    p =".\\reduced\\"+os.path.splitext(os.path.basename(outfile))[0]+".ply"
    print(f"path:{p}")
    pv.PolyData.save(save_cloud,p)

    t1 = o3d.io.read_point_cloud(p)
    t2 = o3d.io.read_point_cloud(p)
    avgd = t1.compute_point_cloud_distance(t2)
    max = np.max(avgd)
    mean = np.mean(avgd)
    print(f"********************************************meax:{max} , mean:{mean}")




#pyvist decimate mesh
def pv_decimate_mesh_NN(infile):
      #las = laspy.read(infile)
    outfile = infile.replace(".ply","_reduced.ply")
    
    print("read done")
    # Create a PyVista PolyData object from the points
    cloud = pv.read(infile)
    scale_factor = 1
    #this is to ensure the plane is created in the correct direction
   # cloud = cloud.scale([scale_factor,scale_factor,1])


   
    tolerance = 0.02
    #cloud.plot(point_size=2)
    cloud2 = cloud.clean(
    point_merging=True,
    merge_tol=tolerance,
    lines_to_points=False,
    polys_to_lines=False,
    strips_to_polys=False,
    inplace=False,
    absolute=False,
    progress_bar=True,  
    )
    #decimated = surf.decimate(0.999)

   


    num_ver = len(cloud2.points)


   



    save_cloud =  pv.PolyData(cloud2.points)
    p =".\\reduced\\"+os.path.splitext(os.path.basename(outfile))[0]+".ply"
    print(f"path:{p}")
    pv.PolyData.save(save_cloud,p)

  

   


#pyvist decimate mesh
def pv_decimate_mesh_visualize(infile):
    #las = laspy.read(infile)
    outfile = infile.replace(".ply","_reduced.ply")
    
    print("read done")
    # Create a PyVista PolyData object from the points
    cloud = pv.read(infile)
    scale_factor = 4
    cloud = cloud.scale([scale_factor,scale_factor,1])

    print(f"number of points {cloud.number_of_points}")
#     tolerance = 0.01
#    # cloud.plot(point_size=2)
#     cloud2 = cloud.clean(
#     point_merging=True,
#     merge_tol=tolerance,
#     lines_to_points=False,
#     polys_to_lines=False,
#     strips_to_polys=False,
#     inplace=False,
#     absolute=False,
#     progress_bar=True,  
#     )

    #print(f"number of points after reduction {cloud2.number_of_points}")

    print("loaded")
    p = pv.Plotter()
    p.add_mesh(cloud,color='red',point_size=1)
    p.show()
    surf = cloud.delaunay_2d(progress_bar=True)

    reduced1 = surf.decimate(0.99,progress_bar=True)
    #decimated = surf.decimate(0.999)

    p = pv.Plotter()
    p.add_mesh(cloud,color='red',point_size=1)
    p.add_mesh(surf)
    p.show()
    p = pv.Plotter()
    p.add_mesh(cloud,color='red',point_size=1)
    p.add_mesh(reduced1)
    p.show()

    num_ver = len(reduced1.points)
    print(f"number of vertexs in the mesh is {num_ver}")

    cloud = cloud.scale([1/scale_factor,1/scale_factor,1])
    reduced1 = reduced1.scale([1/scale_factor,1/scale_factor,1])

    p = pv.Plotter()
    p.add_mesh(cloud,color='red',point_size=1)
    p.add_mesh(reduced1)
    p.show()


    save_cloud =  pv.PolyData(reduced1.points)
    p =".\\reduced\\"+os.path.splitext(os.path.basename(outfile))[0]+".ply"
    print(f"path:{p}")
    pv.PolyData.save(save_cloud,p)

    # surf = cloud.delaunay_2d()
    # decimated = surf.decimate(0.75)

    # wrap3d = pv.PolyData(reduced1.points)
    # wrap3d = pv.wrap(wrap3d)
    # wrap3d.plot()

    # p = pv.Plotter()
    # #p.add_mesh(cloud)
    # p.add_mesh(decimated)
    # p.show()

    



#pyvist decimate mesh
def pv_visualize_mesh(infile):
    #las = laspy.read(infile)
    #outfile = infile.replace(".ply","_reduced.ply")
    
    print("read done")
  
    cloud = pv.read(infile)
    # scale_factor = 4
    # cloud = cloud.scale([scale_factor,scale_factor,1])

    print(f"number of points {cloud.number_of_points}")


    print("loaded")
    p = pv.Plotter()
    p.add_mesh(cloud,color='red',point_size=1)
    p.show()
    surf = cloud.delaunay_2d(progress_bar=True)

 
    p = pv.Plotter()
  #  p.add_mesh(cloud,color='red',point_size=1)
    p.add_mesh(surf)
    p.show()
    p = pv.Plotter()
 


    



#reduce_2_combine("regular3.ply","regular4.ply")
def reduce_2_combine(in_file_1, in_file_2,point_cloud):
    cloud1 = pv.read(in_file_1)
    cloud2 = pv.read(in_file_2)
    print(f"number of points in cloud1 {cloud1.number_of_points}")
    print(f"number of points in cloud2 {cloud2.number_of_points}")
    print("loaded")
    surf1 = cloud1.delaunay_2d(progress_bar=True)

    reduced1 = surf1.decimate(0.99,progress_bar=True)

    surf2 = cloud2.delaunay_2d(progress_bar=True)

    reduced2 = surf2.decimate(0.99,progress_bar=True)

    p = pv.Plotter()
    p.add_mesh(cloud1,color='red',point_size=1)
    p.add_mesh(reduced1)
    p.add_mesh(cloud2,color='yellow',point_size=1)
    p.add_mesh(reduced2)
    p.show()
    combined = pv.PointSet(reduced2.points)
    combined2 = pv.PointSet(reduced1.points)
    combined = combined+combined2

    surf2 = combined.delaunay_2d(progress_bar=True)

    p = pv.Plotter()
   # p.add_mesh(cloud1,color='red',point_size=1)
    p.add_mesh(combined)
    p.add_mesh(surf2)
   # p.add_mesh(cloud2,color='yellow',point_size=1)
    
    p.show()




    

def split_point_cloud(n,infile):
    i = 0
    with laspy.open(infile) as input_las:
        for cur_chunk in input_las.chunk_iterator(10000000):
            point_data = np.stack([cur_chunk.X, cur_chunk.Y, cur_chunk.Z], axis=0).transpose((1, 0))
            point_cloud = o3d.geometry.PointCloud()
            point_cloud.points = o3d.utility.Vector3dVector(point_data)
            fn = os.path.splitext(os.path.basename(infile))[0]
            print(f"./split/{fn}_{i}X.ply")
           # o3d.visualization.draw_geometries([point_cloud]) 
            o3d.io.write_point_cloud(f"./split/{fn}_{i}X.ply", point_cloud)
            i = i + 1

    return



def split_pointcloud_x(n,infile):
    #load the point cloud from a PLY file
   # las = laspy.read(infile)



    point_cloud = o3d.io.read_point_cloud(infile)

    bbox = point_cloud.get_axis_aligned_bounding_box()
    print(f"Bounding box: {bbox}")

    # #calculate dimensions of the bounding box
    x_dim = bbox.max_bound[0] - bbox.min_bound[0]
    y_dim = bbox.max_bound[1] - bbox.min_bound[1]
    z_dim = bbox.max_bound[2] - bbox.min_bound[2]

    #calculate step sizes for splitting
    step_x = x_dim / n
    step_y = y_dim / n

    # #loop through to create n^2 subdivisions
    for i in range(n):
    #     for j in range(n):
            min_bound = np.array([
                bbox.min_bound[0] + i * step_x,
                bbox.min_bound[1],
                bbox.min_bound[2]
            ])
            max_bound = np.array([
                bbox.min_bound[0] + (i + 1) * step_x,
                bbox.max_bound[1],
                bbox.max_bound[2]
            ])
            cropped_pcd = point_cloud.crop(o3d.geometry.AxisAlignedBoundingBox(min_bound, max_bound))
            fn = os.path.splitext(os.path.basename(infile))[0]
            o3d.io.write_point_cloud(f"./split2/{fn}_{i}X.ply", cropped_pcd)
            #pv_visualize_mesh(f"./split2/{fn}_{i}X.ply")
            
    #         #visualize each cropped point cloud
   # o3d.visualization.draw_geometries([cropped_pcd])
            #save each cropped point cloud



def split_dir(in_dir):
    for filename in os.listdir(in_dir):
        f = os.path.join(in_dir, filename)
        # checking if it is a file
        if os.path.isfile(f):
            print(
                f
            )
            split_point_cloud(11,f)
          

def combine_dir(in_dir):
    flist = os.listdir(in_dir)
    f = os.path.join(in_dir, in_dir[1])
    print(f"{f} in { in_dir[0]}")
    curr_point_cloud = o3d.io.read_point_cloud(f)
    for filename in flist[2:]:
        f = os.path.join(in_dir, filename)
        # checking if it is a file
        if os.path.isfile(f):
            print(
                f
            )
            curr_point_cloud = curr_point_cloud+ o3d.io.read_point_cloud(f)
            
        #point_cloud = o3d.io.read_point_cloud(infile)
    o3d.io.write_point_cloud("combined.ply",curr_point_cloud)

def visualize_point_cloud_1():
    i = 0
    with laspy.open("C:\\Users\\code8\\Downloads\\New folder (8)\\split\\RS000011_5X.ply") as input_las:
        for cur_chunk in input_las.chunk_iterator(5000000):
            point_data = np.stack([cur_chunk.X, cur_chunk.Y, cur_chunk.Z], axis=0).transpose((1, 0))
            point_cloud = o3d.geometry.PointCloud()
            point_cloud.points = o3d.utility.Vector3dVector(point_data)
           
          
def split_x_dir(in_dir):
    for filename in os.listdir(in_dir):
        f = os.path.join(in_dir, filename)
        # checking if it is a file
        if os.path.isfile(f):
            print(
                f
            )
            split_pointcloud_x(10,f)

def reduce_dir(in_dir):
     r_sample =  os.listdir(in_dir)
     r_sample = random.sample(r_sample,10)
     for filename in r_sample:
        f = os.path.join(in_dir, filename)
        # checking if it is a file
        if os.path.isfile(f):
            print(
                f
            )
            p =".\\reduced\\"+os.path.splitext(os.path.basename(f))[0]+"_reduced"+".ply"
            if not os.path.exists(p):
                print(f"{p} does not exist")
               # try:
                pv_decimate_mesh(f)
               # except:
                    #print(f"error processing {p}")
            else:
                print(f"{p} already proccesed")

if __name__=="__main__":

    os.chdir("D:\\data_475")

    #split_dir("C:\\Users\\code8\\Downloads\\475_in")
    split_dir("../data2")

    split_x_dir("./split")

    reduce_dir("split2")

    combine_dir("./reduced")
   
    #visualizes the data
   # visualize_point_cloud_1()
    pv_visualize_mesh("./combined.ply")
  
