# Generate a PLSG with L-shape with random points inside of
# Input: nPoints: number of points inside the square and in the borders
#        Tolerance: Distance to consirate a point to close to the border, so it must move to the border
# Output: .poly file


import numpy
import random
import sys
import math
from shapely.geometry import Point, Polygon, LineString
from shapely.ops import nearest_points
import sys

def insert_border_ccw(border, pt):
    for i in range(len(border)):
        e1 = border[i % len(border)]
        e2 = border[(i+1) % len(border)]
        line = LineString([Point(e1),Point(e2)])
        #print(e1,e2,line.wkt)
        if (pt.within(line)):
            border.insert((i+1)%len(border),(pt.x, pt.y))
            break

def square2x2(num_points, tolerance, percentage_center_point = 0.1, radius_center = 0.2):
    border_original = [(-1.0, 1.0), (-1.0, -1.0), (0.0, -1.0), (0.0, 0.0), (1.0, 0.0), (1.0,1.0)]
    new_border = border_original.copy()

    poly = Polygon(border_original) 
    
    s = set()

    points_in_center = int(num_points*percentage_center_point)
    #insert in center
    while len(s) + len(new_border) - 5 != points_in_center:
        x = random.uniform(-1.0*radius_center, 1.0*radius_center )
        y = random.uniform(-1.0*radius_center, 1.0*radius_center )
        pt = Point(x,y)
        
        if(poly.contains(pt) and numpy.sqrt(x*x + y*y) <= (radius_center)*(radius_center) ):
            #distance = numpy.round(poly.boundary.distance(pt), len(str(tolerance))-1 )
            if( poly.boundary.distance(pt) < tolerance*radius_center):
                border_pt = poly.exterior.interpolate(poly.boundary.project(pt))
                #border_pt = nearest_points(poly, pt)
                #s.add((border_pt[0].x, border_pt[0].y))
                insert_border_ccw(new_border, border_pt)
            else:
                s.add((x, y))

    insert_border_ccw(new_border, Point(0.02, 0.2))
    insert_border_ccw(new_border, Point(0.02, 0.0))
    insert_border_ccw(new_border, Point(0.01, 0.0))
    insert_border_ccw(new_border, Point(0.001, 0.0))
    insert_border_ccw(new_border, Point(0.0001, 0.0))

    

    while len(s) + len(new_border) != num_points:
        x = random.uniform(-1.0, 1.0 )
        y = random.uniform(-1.0, 1.0 )
        pt = Point(x,y)
        if(poly.contains(pt) and numpy.sqrt(x*x + y*y) > (radius_center)*(radius_center) ):
            #distance = numpy.round(poly.boundary.distance(pt), len(str(tolerance))-1 )
            if( poly.boundary.distance(pt) < tolerance):
                border_pt = poly.exterior.interpolate(poly.boundary.project(pt))
                #border_pt = nearest_points(poly, pt)
                #s.add((border_pt[0].x, border_pt[0].y))
                insert_border_ccw(new_border, border_pt)
            else:
                s.add((x, y))

    write_node(new_border + list(s),tolerance, percentage_center_point, radius_center)



def write_node(points, tolerance, percentage_center_point, radius_center):
    n_points =len(points)
    print("Printing .node file with {} points in {}".format(n_points, "LR" + str(n_points) + '.node'))
    f = open("LR" + str(n_points) + '.node', 'w')
    f.write("# Points: {}, tolerance: {}, Percentage center point: {}, radius_center: {}\n".format(n_points, tolerance, percentage_center_point, radius_center))
    f.write("{} 2 0 0\n".format(n_points))
    for i in range(0, n_points):
        f.write('{0} {1} {2}\n'.format(i, points[i][0], points[i][1]))
    f.close()

if __name__ == "__main__":
    random.seed(138)
    full_cmd_arguments = sys.argv
    argument_list = full_cmd_arguments[1:]
    num_points = int(argument_list[0])
    tolerance = tolerance = 1/10**(0.5*math.log10(num_points) + 0.5)
    square2x2(num_points, tolerance)
