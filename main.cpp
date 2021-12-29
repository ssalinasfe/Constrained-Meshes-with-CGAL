#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/IO/Color.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Constrained_triangulation_2.h>
#include <CGAL/draw_triangulation_2.h>
#include <fstream>

#include <iterator>
#include <array>
#include <vector>
#include <iostream>
#include <cmath>
#include <string>
#include <sstream>
#include <random>       // std::default_random_engine
#include <algorithm>    // std::shuffle
#include <cstdlib> 
#include <cctype>
#include <iomanip>

struct ElementInfo
{
  int index;
  bool border = false;
};


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_with_info_2<ElementInfo, K>    Vb;
typedef CGAL::Triangulation_face_base_with_info_2<ElementInfo,K> Fbb;//era fb
typedef CGAL::Constrained_triangulation_face_base_2<K,Fbb> Fb;  //added
typedef CGAL::Triangulation_data_structure_2<Vb,Fb> Tds;
typedef CGAL::Exact_predicates_tag Itag; // added
typedef CGAL::Constrained_triangulation_2<K, Tds, Itag>  Triangulation;
typedef Triangulation::Constrained_edges_iterator Constrained_edges_iterator;
typedef Triangulation::Face_handle Face_handle;
typedef Triangulation::Edge Edge;
typedef Triangulation::Point  Point;
typedef CGAL::Polygon_2<K> Polygon_2;

//Read node file in .node format and nodes in point vector
std::vector<Point> read_nodes_from_file(std::string name){
    std::string line;
    std::ifstream nodefile(name);
    double a1, a2, a3, a4;
    std::vector<Point> points;
    unsigned n_vertices = 0;
    std::string tmp;
    if (nodefile.is_open())
    {
      while (std::getline(nodefile, line)){ //add check boundary vertices flag
        std::istringstream(line) >> tmp;
        if (tmp[0] != '#' ) //check if first line is a comentary
        {
            n_vertices = std::stoi(tmp);
            points.reserve(n_vertices);
            break;
        }
      }
      while (std::getline(nodefile, line))
      {
        std::istringstream(line) >> tmp;
        if (tmp[0] != '#' ) //check if first line is a comentary
        {
          std::istringstream(line) >> a1 >> a2 >> a3 >> a4;
          points.push_back(Point(a2, a3));
        }
      }  
    }
    else 
        std::cout << "Unable to open node file"; 
    nodefile.close();
    return points;
}


void print_off(Triangulation Tr, std::string name, int n_faces){
    std::ofstream outfile(name + ".1.off");
    outfile << "OFF" << std::endl;
    outfile << Tr.number_of_vertices() << " " << n_faces << " " << 0 << std::endl;
    for(Triangulation::Vertex_handle v : Tr.finite_vertex_handles()){
        outfile << std::setprecision(15)<< v->point().x() <<" "<<v->point().y()<<" 0.0\n";
    }
    for(Triangulation::Face_handle f : Tr.finite_face_handles()){
      if(f->info().index!=-1)
        outfile << 3 << " " << f->vertex(0)->info().index << " " << f->vertex(1)->info().index << " " << f->vertex(2)->info().index << std::endl;
    }
    outfile.close();
}


void print_node_file(Triangulation Tr, std::string name){
    std::ofstream outfile(name + ".1.node");
    outfile << Tr.number_of_vertices() <<" 2 0 1"<< std::endl;
    for(Triangulation::Vertex_handle v : Tr.finite_vertex_handles()){
        outfile << v->info().index<<" "<< std::setprecision(15)<< v->point().x() << " " << v->point().y() << " "<<v->info().border << std::endl;
    }
    outfile.close();
}

void print_ele_file(Triangulation Tr, std::string name, int n_faces){
    std::ofstream outfile(name + ".1.ele");
    outfile << n_faces <<" 3 0"<< std::endl;
    for(Triangulation::Face_handle f : Tr.finite_face_handles()){
        if(f->info().index!=-1)
            outfile << f->info().index << " " << f->vertex(0)->info().index << " " << f->vertex(1)->info().index << " " << f->vertex(2)->info().index << std::endl;
    }
    outfile.close();
}

void print_neigh_file(Triangulation Tr, std::string name, int n_faces){
  std::ofstream outfile(name + ".1.neigh");
  outfile << n_faces <<" 3"<< std::endl;
  for(Triangulation::Face_handle f : Tr.finite_face_handles()){
      if(f->info().index!=-1){
          outfile << f->info().index << " " << f->neighbor(0)->info().index << " " << f->neighbor(1)->info().index << " " << f->neighbor(2)->info().index << std::endl;
      }
  }
  outfile.close();
}

void print_neigh_file(Triangulation Tr){
  for(Triangulation::Finite_faces_iterator fit = Tr.finite_faces_begin(); fit != Tr.finite_faces_end(); ++fit)
    if(fit->info().index != -1)
      std::cout << "neigh " << fit->info().index << ": (" << fit->neighbor(0)->info().index << ", " << fit->neighbor(1)->info().index << ", " << fit->neighbor(2)->info().index << ")" << std::endl;
}


bool check_inside(Point pt, std::vector<Point> boundary, K traits)
{
  switch(CGAL::bounded_side_2(boundary.begin(), boundary.end(), pt, traits)) {
    case CGAL::ON_BOUNDED_SIDE :
      return true;
    case CGAL::ON_BOUNDARY:
      return true;
    case CGAL::ON_UNBOUNDED_SIDE:
      return false;
  }
  return false;
}



void printTriangulation(Triangulation Tr, std::vector<Point> boundary, std::string output_file){
  int index = 0;
  for(Triangulation::Vertex_handle v : Tr.finite_vertex_handles()) {
    v->info().index = index;
    index++;
  }

  for (Constrained_edges_iterator it = Tr.constrained_edges_begin(); it != Tr.constrained_edges_end(); ++it) {
    Edge e = *it;
    e.first->vertex( (e.second+1)%3 )->info().border=true;
    e.first->vertex( (e.second+2)%3 )->info().border=true;
  }

  index = 0;
  for(Triangulation::Face_handle f : Tr.all_face_handles()){
    if(!Tr.is_infinite(f)){
      auto v0 = f->vertex(0)->point();
      auto v1 = f->vertex(1)->point();
      auto v2 = f->vertex(2)->point();
      Point mid1 = Point((v0.x()+v1.x())/2,(v0.y()+v1.y())/2);
      Point mid2 = Point((v1.x()+v2.x())/2,(v1.y()+v2.y())/2);
      Point mid3 = Point((v2.x()+v0.x())/2,(v2.y()+v0.y())/2);
      if(check_inside(mid1, boundary, K()) && check_inside(mid2,boundary, K()) && check_inside(mid3, boundary, K())){
        f->info().index = index;
        index++;
      }
      else{
        f->info().index = -1;
      }
    }else{
      f->info().index = -1;
    }
  }

  std::cout<<"Printing mesh with "<<index<<" faces and "<< Tr.number_of_vertices() <<" vertices"<<std::endl;
  print_off(Tr, output_file, index);
  print_node_file(Tr, output_file);
  print_ele_file(Tr, output_file, index);
  print_neigh_file(Tr, output_file, index);
  std::cout<<"Printed files in "<<output_file<<std::endl;
}

int main(int argc, char **argv) {

  //Read nodes and output file name
  std::vector<Point> points = read_nodes_from_file(argv[1]); 
  std::string output_file = std::string(argv[2]);


  std::cout<<"Read "<<points.size()<<" points from file "<<argv[1]<<std::endl;
  unsigned seed = 138;

  //Shuffle points
  std::shuffle( points.begin(), points.end(), std::default_random_engine(seed));
  std::vector<Point> boundary = { Point(0.0,0.0), Point(1.0,0.0), Point(1.0,1.0), Point(-1.0,1.0), Point(-1.0,-1.0), Point(0.0,-1.0) };
  


  //Triangulation generation
  Triangulation Tr;
  Tr.insert(points.begin(),points.end());
  Tr.insert_constraint(boundary.begin(), boundary.end(), true);
  printTriangulation(Tr, boundary, output_file);
  

  //CGAL::draw(Tr);
  return EXIT_SUCCESS;
}