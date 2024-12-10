import laspy
import open3d as o3d
import numpy as np

import pyvista as pv
#yes
import os as os
import random




def split_point_cloud_n(n,inf):

  

    point_cloud_in =  o3d.io.read_point_cloud(inf)

    point_cloud  = o3d.geometry.PointCloud()
    point_cloud.points = o3d.utility.Vector3dVector(point_cloud_in.points)



    bbox = point_cloud.get_axis_aligned_bounding_box()
    print(f"Bounding box: {bbox}")

    #calculate dimensions of the bounding box
    x_dim = bbox.max_bound[0] - bbox.min_bound[0]
    y_dim = bbox.max_bound[1] - bbox.min_bound[1]
    z_dim = bbox.max_bound[2] - bbox.min_bound[2]

    #calculate step sizes for splitting
    step_x = x_dim / n
    step_y = y_dim / n

    #loop through to create n^2 subdivisions
    for i in range(n):
        for j in range(n):
            min_bound = np.array([
                bbox.min_bound[0] + i * step_x,
                bbox.min_bound[1] ,
                bbox.min_bound[2]
            ])
            max_bound = np.array([
                bbox.min_bound[0] + (i + 1) * step_x,
                bbox.max_bound[1],
                bbox.max_bound[2]
            ])
            cropped_pcd = point_cloud.crop(o3d.geometry.AxisAlignedBoundingBox(min_bound, max_bound))
            
           
            #save each cropped point cloud
            print(f"cropped_{i}X{j}.ply")
            o3d.io.write_point_cloud(f"cropped_{i}X{j}.ply", cropped_pcd)
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
   # p.show()
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


    surf = cloud2.delaunay_2d(progress_bar=True)

    #reduced1 = surf.decimate(0.95,progress_bar=True)
   
    p = pv.Plotter()
    p.add_mesh(cloud,color='red',point_size=1)
    p.add_mesh(surf)
    p.show()

    num_ver = len(cloud2.points)


   



    save_cloud =  pv.PolyData(cloud2.points)
    p =".\\reduced\\"+os.path.splitext(os.path.basename(outfile))[0]+".ply"
    print(f"path:{p}")
    pv.PolyData.save(save_cloud,p)
    

def pv_reduce_and_compare(infile):
    outfile = infile.replace(".ply","_reduced.ply")
    
    print("read done")
    # Create a PyVista PolyData object from the points
    cloud = pv.read(infile)
    scale_factor = 4
    #this is to ensure the plane is created in the correct direction
    cloud = cloud.scale([scale_factor,scale_factor,1])


   
    tolerance = 0.01
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

    ratio = cloud2.number_of_points/cloud.number_of_points
    #decimated = surf.decimate(0.999)


    

    surf_cloud = cloud2.delaunay_2d(progress_bar=True)
    surf_cloud = surf_cloud.scale([1/scale_factor,1/scale_factor,1])
    surf2 = cloud.delaunay_2d(progress_bar=True)
    cloud = cloud.scale([1/scale_factor,1/scale_factor,1])
    surf2 = surf2.scale([1/scale_factor,1/scale_factor,1])
    print(f"ratio: {1-ratio}")
    surf3 = surf2.decimate(1-ratio,progress_bar=True)
    

    p = pv.Plotter(shape=(1,2))
  
   
    p.subplot(0, 0)
    p.add_text("cloud reduction", font_size=30)
    p.add_mesh(surf_cloud)
    p.add_mesh(cloud,color="red", point_size=1)
    p.subplot(0, 1)
    p.add_text("mesh reduction", font_size=30)
    p.add_mesh(surf3)
    p.add_mesh(cloud,color="red", point_size=1)
    p.link_views()
    p.show()


    num_ver = len(cloud2.points)

    


   



    save_cloud =  pv.PolyData(cloud2.points)
    p =".\\reduced\\"+os.path.splitext(os.path.basename(outfile))[0]+".ply"
    print(f"path:{p}")

    #surf4 = cloud.reconstruct_surface(progress_bar=True)


    # p = pv.Plotter()
    
    # slices = surf_cloud.slice_along_axis(n=10, axis='z')
    
    # slices2 = surf3.slice_along_axis(n=10, axis='z')
    # p.add_mesh(slices, color="green")
    # p.add_mesh(slices2, color="blue")
    # #p.add_mesh(surf4)
    # p.show()

    # p = pv.Plotter()
    # p.add_mesh(surf_cloud)
    # p.add_mesh(surf3)
    # p.show()
    
    pv.PolyData.save(save_cloud,p)



def split_dir(in_dir):
    for filename in os.listdir(in_dir):
        f = os.path.join(in_dir, filename)
        # checking if it is a file
        if os.path.isfile(f):
            print(
                f
            )
            split_point_cloud(11,f)


os.chdir("D:\\data_475\\")
l = os.listdir("./split2")

l = list(l)
l = random.sample(l,5)

for filename in l:
    f = os.path.join("./split2", filename)
    # checking if it is a file
    if os.path.isfile(f):
        print(
            f
        )
    
        pv_reduce_and_compare(f)

