#include <iostream>
#include <opencv2/opencv.hpp>
#include <list>

std::map<int, cv::Vec3b> colors;

namespace Geometry {

	typedef cv::Point2d Point;
	typedef std::vector<Point> Polygon;
	typedef std::array<Point,3> Triangle;

	Polygon cut_polygon;


	// Given three collinear points p, q, r, the function checks if
	// point q lies on line segment 'pr'
	bool on_segment(Point p, Point q, Point r) {
		if (q.x < std::max(p.x, r.x) && q.x > std::min(p.x, r.x) &&
			q.y < std::max(p.y, r.y) && q.y > std::min(p.y, r.y))
		   return true;
	 
		return false;
	}
	 
	// To find orientation of ordered triplet (p, q, r).
	// The function returns following values
	//  0 --> p, q and r are collinear
	//  1 --> Clockwise
	// -1 --> Counterclockwise
	int orientation(Point p, Point q, Point r) {
		double v = (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y);
		if (v>0) return  1;
		if (v<0) return -1;
		return 0;
	}
	 
	// The main function that returns true if line segment 'p1q1'
	// and 'p2q2' intersect.
	bool do_intersect(Point p1, Point q1, Point p2, Point q2) {
		// Find the four orientations needed for general and
		// special cases
		double o1 = orientation(p1, q1, p2);
		double o2 = orientation(p1, q1, q2);
		double o3 = orientation(p2, q2, p1);
		double o4 = orientation(p2, q2, q1);
	 
		// General case
		if (o1 != o2 && o3 != o4)
			return true;
	 
		// Special Cases
		// p1, q1 and p2 are collinear and p2 lies on segment p1q1
		if (o1 == 0 && on_segment(p1, p2, q1)) return true;
	 
		// p1, q1 and q2 are collinear and q2 lies on segment p1q1
		if (o2 == 0 && on_segment(p1, q2, q1)) return true;
	 
		// p2, q2 and p1 are collinear and p1 lies on segment p2q2
		if (o3 == 0 && on_segment(p2, p1, q2)) return true;
	 
		 // p2, q2 and q1 are collinear and q1 lies on segment p2q2
		if (o4 == 0 && on_segment(p2, q1, q2)) return true;
	 
		return false; // Doesn't fall in any of the above cases
	}


	std::vector<Point> make_polygon_clockwise(Polygon polygon) {
		
		int pt = 0;
		for (uint i=1; i<polygon.size(); i++) {
			
			if (polygon[i].x < polygon[pt].x) {
				pt = i;
			} else if (polygon[i].x == polygon[pt].x && polygon[i].y < polygon[pt].y) {
				pt = i;
			}
		}
		
		pt = pt + polygon.size();
		
		auto a = polygon[(pt-1)%polygon.size()];
		auto b = polygon[(pt  )%polygon.size()];
		auto c = polygon[(pt+1)%polygon.size()];
		
		if ((b.x*c.y+a.x*b.y+a.y*c.x)-(a.y*b.x+b.y*c.x+a.x*c.y) < 0) {
			
			std::reverse(polygon.begin(), polygon.end());
			std::cerr << "Reversing polygon" << std::endl;
		}

		return polygon;
	}

	bool is_point_in_triangle(Point p, Triangle tr) {
		
		Point p0 = tr[0], p1 = tr[1], p2 = tr[2];
		
		double s = (p0.x - p2.x) * (p.y - p2.y) - (p0.y - p2.y) * (p.x - p2.x);
		double t = (p1.x - p0.x) * (p.y - p0.y) - (p1.y - p0.y) * (p.x - p0.x);

		if ((s < 0) != (t < 0) && s != 0 && t != 0)
			return false;

		double d = (p2.x - p1.x) * (p.y - p1.y) - (p2.y - p1.y) * (p.x - p1.x);
		return d == 0 || (d < 0) == (s + t <= 0);
	}

	bool is_point_in_triangle_circle(Point p, Triangle tr) {
		
		Point p0 = tr[0], p1 = tr[1], p2 = tr[2];
		
		if ((p1.x*p2.y+p0.x*p1.y+p0.y*p2.x)-(p0.y*p1.x+p1.y*p2.x+p0.x*p2.y) < 0) {
			
			//std::swap(p0,p1);
			std::cerr << "Reversing triangle" << std::endl;
		}

		
		double m00 =  p0.x-p.x, m01 = p0.y-p.y, m02 = p0.x*p0.x-p.x*p.x+p0.y*p0.y-p.y*p.y;
		double m10 =  p1.x-p.x, m11 = p1.y-p.y, m12 = p1.x*p1.x-p.x*p.x+p1.y*p1.y-p.y*p.y;
		double m20 =  p2.x-p.x, m21 = p2.y-p.y, m22 = p2.x*p2.x-p.x*p.x+p2.y*p2.y-p.y*p.y;

		double d = m00*m11*m22 + 
				m01*m12*m20 + 
				m02*m10*m21 - 
				m20*m11*m02 - 
				m10*m01*m22 - 
				m00*m21*m12;
		
		return d > 0;
	}
	
	std::vector<Triangle> delaunay(std::vector<Triangle> triangles) {
	
		bool found = true;
		while (found) {
			found = false;
			for (uint i=0; i<triangles.size(); i++) {
				for (uint j=0; j<i; j++) {
					
					bool found2 = false;
					for (uint ii=0; not found2 and ii<3; ii++) {
						for (uint jj=0; not found2 and jj<3; jj++) {
							if (triangles[i][0].x == triangles[j][1].x and
								triangles[i][0].y == triangles[j][1].y and
								triangles[i][1].x == triangles[j][0].x and
								triangles[i][1].y == triangles[j][0].y) {
								
								found2 = true;
								break;
							}
							
							std::swap( triangles[i][0], triangles[i][1]); 
							std::swap( triangles[i][0], triangles[i][2]); 

						}
						if (not found2) {
							std::swap( triangles[j][0], triangles[j][1]); 
							std::swap( triangles[j][0], triangles[j][2]); 
						}
					}
					
					if (found2) {

						if ( 
							is_point_in_triangle_circle( triangles[j][2], triangles[i] ) or
							is_point_in_triangle_circle( triangles[i][2], triangles[j] ) ) {
						
							Point a = triangles[i][0], b = triangles[i][1], c = triangles[i][2], d = triangles[j][2];
							
							triangles[i][0] = a;
							triangles[i][1] = d;
							triangles[i][2] = c;

							triangles[j][0] = d;
							triangles[j][1] = b;
							triangles[j][2] = c;
							
							found = true;
						}	
					}				
				}
			}
		}
		
		return triangles;		
	}

	std::vector<Triangle> triangulate(Polygon polygon) {
		
		polygon = make_polygon_clockwise(polygon);
		
		std::vector<Triangle> triangles;
		{
			std::list<Point> mp;
			for (auto &p : polygon) mp.push_back(p);
			
			while (mp.size()>=3) {
				
				{
					auto ita = mp.begin(), itb = mp.begin(), itc = mp.begin();
					++itb;
					++++itc;			
					bool found = false;
					while (ita!=mp.end()) {
						
						Triangle t = { *ita, *itb, *itc };
						
						if (mp.size()==3) {

							mp.clear();
							found = true;
							triangles.push_back(t);
							break;				
						}
							
						if ((t[1].x*t[2].y+t[0].x*t[1].y+t[0].y*t[2].x)-(t[0].y*t[1].x+t[1].y*t[2].x+t[0].x*t[2].y) > 0) {
							
							auto itd = itc;
							
							++itd;
							
							bool found2 = false;
							while (itd!=ita and not found2) {
								
								auto d = *itd;
								
								if (d != t[0] and d != t[1] and d != t[2] ) {
									found2 = is_point_in_triangle(d,t);
								}
								
								
								++itd;
								if (itd==mp.end()) itd=mp.begin();
							}
							
							//std::cerr << "_*"[found2];
							
							if (not found2) {
							
								mp.erase(itb);
								found = true;
								triangles.push_back(t);
								break;
							}
						}
						
						ita++;
						itb++;
						itc++;
						if (itb==mp.end()) itb=mp.begin();
						if (itc==mp.end()) itc=mp.begin();
					}
					if (not found) {
						
						std::cerr << "Triangulation not found" << std::endl;
						return triangles;
					}
				}
			}
		}
		return triangles;
	}
	
	std::vector<Triangle> triangulate_polygon_with_holes(Polygon polygon, std::vector<Polygon> holes) {
		
		
		if (polygon.size()<3) return std::vector<Triangle>();
		
		std::list<Point> mp;
		mp.push_back(polygon.back());
		for (auto &p : polygon) {
			
			double distance = cv::norm((mp.back()-p));
			int steps = distance / 64;
			auto step_vector = (p-mp.back())/(1+steps);
			
			for (int i=0; i<steps; i++) {
				
				auto noise = [](){return (rand()%(1<<12) - (1<<11)) / double(1<<26); };
				
				mp.push_back(mp.back() + step_vector + Point(noise(), noise()));
			}
			
			mp.push_back(p);
		}
		mp.pop_front();
		
		for (auto &hole : holes) {
			
			if (hole.size()<3) continue;

			hole = make_polygon_clockwise(hole);
			std::reverse(hole.begin(), hole.end());
			
			bool found_insert = false;
			for (auto ith = hole.begin(); not found_insert and ith != hole.end(); ith++) {
				
				auto &ph = *ith;
				
				for (auto it = mp.begin(); not found_insert and it!=mp.end(); it++) {
					
					auto a = ph, b = *it;
					bool found_intersection = false;
					
					for (auto itr = mp.begin(); not found_intersection and itr!=mp.end(); itr++) {
						
						found_intersection = itr != it and *itr == *it;
					}
					
					for (auto it2 = mp.begin(); not found_intersection and it2!=mp.end(); it2++) {
						
						auto it3 = it2;
						it3++;
						if (it3 == mp.end()) it3 = mp.begin();
						
						if (b==*it2) continue;
						if (b==*it3) continue;
						
						found_intersection = do_intersect(a, b, *it2, *it3);
					}
					
					for (auto &hh : holes) {
						
						for (int i = hh.size()-1; not found_intersection and i>=0; i--) {
							auto c = hh[i];
							auto d = hh[(i+1)%hh.size()];
							
							
							if (a==c or a==d or b==c or b==d) continue;

							found_intersection = do_intersect(a, b, c, d);
						}
					}
					
					if (not found_intersection) {
						
						
						mp.insert(it, *it);

						for (auto itj = ith; itj < hole.end(); itj++) 
							mp.insert(it, *itj);
						for (auto itj = hole.begin(); itj < ith; itj++) 
							mp.insert(it, *itj);

						mp.insert(it, *ith);
						
						found_insert = true;
						//std::cerr << "Inserted!" << std::endl;
						
					}
				}
			}
		}
		
		polygon.clear();
		for (auto &p : mp) polygon.push_back(p);

		std::vector<Triangle> triangles = delaunay(triangulate(polygon));
		//std::vector<Triangle> triangles = triangulate(polygon);
		
		cut_polygon = polygon;
		
		return triangles;
	}

	struct Scenario {
		
		std::vector<Triangle> walkable_areas;
		
		std::vector<std::pair<Triangle, int>> attractors;
	};
	
	Scenario build_scenario(std::vector<std::pair<Polygon, std::vector<Polygon>>> polygons_and_holes, cv::Size sz) {
	
		Scenario scenario;
	
		for (auto &[polygon, holes] : polygons_and_holes) {
			
			auto triangles = triangulate_polygon_with_holes(polygon, holes);
			
			for (auto &t : triangles) {
				
				scenario.walkable_areas.push_back(t);
			}
		}
		
		double scale = 4;
		
		cv::Mat1i colored_field(sz.height*scale, sz.width*scale);
		
		std::map<Point, int> already_colored;
		
		for (int i = 0; i<colored_field.rows; i++) {
			for (int j = 0; j<colored_field.cols; j++) {
				
				colored_field(i,j) = -1;
				
				double max_distance = 1e10;
				int idx = 0;
				for (auto &t : scenario.walkable_areas) {
					
					Point pt = {j/scale,i/scale};
					
					if (is_point_in_triangle(pt, t)) {
						
						colored_field(i,j) = -1;
						break;
					}
					
					for (int ti=0; ti<3; ti++) {
						{
							double d = cv::norm(pt - t[ti]);
							if (d<max_distance) { colored_field(i,j) = idx; max_distance = d; }
							idx++;
						}

						Point p0 = t[ti];
						Point p1 = t[(ti+1)%3];
						
						if ((pt-p0).dot(p1-p0) >= 0) {
							if ((pt-p1).dot(p0-p1) >= 0) {
								double d = ((p1.x-p0.x)*(p0.y-pt.y)-(p0.x-pt.x)*(p1.y-p0.y))/cv::norm(p1-p0);
								if (d>0 and d<=max_distance) { colored_field(i,j) = idx; max_distance = d; }
							}
						}
						idx++;
					}
				}
			}
		}

		if (scenario.walkable_areas.size()) {
			
			cv::Mat3b screen(colored_field.size());
			
			colors[-1] = cv::Vec3b(0,0,0);
			
			for (int i = 0; i<colored_field.rows; i++) {
				for (int j = 0; j<colored_field.cols; j++) {
				
					colors.insert({colored_field(i,j), cv::Vec3b(rand()%256, rand()%256, rand()%256)});
					screen(i,j) = colors[colored_field(i,j)];
				}
			}
			
			
			
			//cv::Mat3b screen44;
			//cv::resize(screen, screen44, cv::Size(), 4, 4, cv::INTER_NEAREST);
			cv::imshow("poligon_debug", screen);
			
		}

		std::map<int, std::vector<cv::Point>> polygons;

		for (int i = 0; i<colored_field.rows; i++) {
			for (int j = 0; j<colored_field.cols; j++) {
				
				cv::Point pt = {j,i};
				
				if (i && j && i!=colored_field.rows-1 && j!=colored_field.cols-1) {
					
					if (colored_field(i,j) != colored_field(i+1,j  ) ||
						colored_field(i,j) != colored_field(i-1,j  ) ||
						colored_field(i,j) != colored_field(i  ,j+1) ||
						colored_field(i,j) != colored_field(i  ,j-1)) {
					
						polygons[colored_field(pt)].push_back(pt);
					}
					
				} else {
					polygons[colored_field(pt)].push_back(pt);
				}
			}
		}
		
		polygons.erase(-1);
		
		std::map<int, std::vector<cv::Point>> hulls;
		for (auto &ph : polygons) {
			std::vector<cv::Point> hull;
			cv::convexHull(ph.second, hull, false);

			std::vector<cv::Point> hull2;	
			for (int i=0; i<hull.size(); i++) {
					
				Point a = hull[i] - hull[(i+1)%hull.size()];
				Point b = hull[(i+1)%hull.size()] - hull[(i+2)%hull.size()];
				
				if (std::abs(a.dot(b)/cv::norm(a)/cv::norm(b)) < 0.95) {
					hull2.push_back(hull[(i+1)%hull.size()]);
				}
			}
			
			for (auto &t : scenario.walkable_areas) {
				
				for (int i=0; i<3; i++) {
					
					for (auto &h : hull2) {
						
						if (cv::norm(t[i]-Point(h)) < 4*scale) 
							h = t[i];
							
						if (h.x < scale) h.x = 0;
						if (h.y < scale) h.y = 0;

						if (h.x > colored_field.cols-scale-1) h.x = colored_field.cols-1;
						if (h.y > colored_field.rows-scale-1) h.y = colored_field.rows-1;
					}
				}
				
			}

			hull = hull2;
			hull2.clear();
			for (int i=0; i<hull.size(); i++) {
					
				Point a = hull[i] - hull[(i+1)%hull.size()];
				
				if (a.dot(a) > scale*scale) {
					hull2.push_back(hull[(i+1)%hull.size()]);
				}
			}
			
			hulls[ph.first] = hull2;
		}
		
		std::vector<Point> hull_cluster_points;
		
		
		
		std::map<std::pair<int, int>, std::set<int>> hull_point_count;
		for (auto &ph : hulls) {
			for (auto &p : ph.second) {
				
				if (p.x>0 && p.x<colored_field.cols-1 && p.y>0 && p.y<colored_field.rows-1) {
					
					std::set<int> colors;
					colors.insert(colored_field(p.y-1,p.x-1));
					colors.insert(colored_field(p.y-1,p.x  ));
					colors.insert(colored_field(p.y-1,p.x+1));

					colors.insert(colored_field(p.y  ,p.x-1));
					colors.insert(colored_field(p.y  ,p.x  ));
					colors.insert(colored_field(p.y  ,p.x+1));

					colors.insert(colored_field(p.y+1,p.x-1));
					colors.insert(colored_field(p.y+1,p.x  ));
					colors.insert(colored_field(p.y+1,p.x+1));
					
					if (colors.size()<3) continue;
				}
				
				hull_point_count[{p.x-1, p.y-1}].insert(ph.first);
				hull_point_count[{p.x-1, p.y  }].insert(ph.first);
				hull_point_count[{p.x-1, p.y+1}].insert(ph.first);

				hull_point_count[{p.x  , p.y-1}].insert(ph.first);
				hull_point_count[{p.x  , p.y  }].insert(ph.first);
				hull_point_count[{p.x  , p.y+1}].insert(ph.first);

				hull_point_count[{p.x+1, p.y-1}].insert(ph.first);
				hull_point_count[{p.x+1, p.y  }].insert(ph.first);
				hull_point_count[{p.x+1, p.y+1}].insert(ph.first);
			}
		}
			
		
		for (auto &ph : hulls) {
			
			if (ph.first==-1) continue;
			
			auto &hull = ph.second;
			Polygon polygon;
			if (false) {
				polygon.push_back(hull.front());
				for (int i=1; i<hull.size()-1; i++) {
					
					Point a = polygon.back() - Point(hull[i]);
					Point b = hull[i+1] - hull[i];
					
					if (std::abs(a.dot(b)/cv::norm(a)/cv::norm(b)) < 0.95) {
						polygon.push_back(hull[i]);
					}
				}
				polygon.push_back(hull.back());
				
				for (auto &p : polygon) p = p / scale;
			}

			for (int i=0; i<hull.size(); i++)
			//	if (hull_point_count[{hull[i].x, hull[i].y}].size() > 1)
					polygon.push_back(hull[i]/scale);
				
			if (polygon.size()<3) continue;
			
			while (polygon.size()<3) polygon.push_back(polygon.back());

			std::cerr << hull.size() << " " << polygon.size() << std::endl;
			
			for (auto t : delaunay(triangulate(polygon))) {
				scenario.attractors.push_back({t, ph.first});
			}
			
			

		}
		
		return scenario;
		
	}

};




Geometry::Polygon main_polygon;
std::vector<Geometry::Polygon> holes;

int status = 0;

void MouseCallback(
	int event,
	int x,
	int y,
	int flags,
	void *
	) {
	
	x /= 4;
	y /= 4;
	
	auto noise = [](){return (rand()%(1<<20) - (1<<19)) / double(1<<26); };
		
	if (event == cv::EVENT_LBUTTONDOWN) {
		
		std::cerr << "EVENT_LBUTTONDOWN " << x << " " << y << std::endl;
		
		if (status==0) {
			
			main_polygon.push_back({x + noise() + 0.25, y + noise() + 0.25});
			
			std::cerr << main_polygon.back().x << " " << main_polygon.back().y << std::endl;
			
		}

		if (status==1) {
			
			holes.back().push_back({x + noise() + 0.25 , y + noise() + 0.25});
		}
		
	}

	if (event == cv::EVENT_RBUTTONDOWN) {
		
		std::cerr << "EVENT_RBUTTONDOWN " << x << " " << y << std::endl;
	}
	
	
}

int main() {

	
	{
		cv::Mat3b screen(192,256);
		
		cv::Mat3b screen44;
		cv::resize(screen, screen44, cv::Size(), 4, 4, cv::INTER_NEAREST);
		cv::imshow("poligon", screen44);
		
		cv::setMouseCallback( "poligon", MouseCallback, nullptr );
	}
    
    while (true) {
		int k = cv::waitKey(100);
		
		

		if (k==32 && status==1) {
			
			if (holes.back().empty()) {
				
				holes.pop_back();
				status = 2;
			
			} else {
				
				holes.back() = Geometry::make_polygon_clockwise(holes.back());
				holes.emplace_back();
			}
		}
		
		if (k==32 && status==0) {
			
				status = 1;
				main_polygon = Geometry::make_polygon_clockwise(main_polygon);
				holes.emplace_back();
		}
		
		cv::Mat3b screen(192,256,cv::Vec3b(0,0,0));
		

		{
			
			Geometry::Scenario scenario = Geometry::build_scenario({{main_polygon, holes}},{256,192});
			auto triangles = scenario.walkable_areas;
			
			for (auto triangle : triangles) {
				for (uint i=0; i<triangle.size(); i++) {
					
					Geometry::Point &ca = triangle[i];
					Geometry::Point &cb = triangle[(i+1) % triangle.size()];
					cv::line(screen, ca, cb, cv::Scalar(64,64,64));
				}			
			}
			
			
			for (auto triangle_pair : scenario.attractors) {
				
				auto triangle = triangle_pair.first;
				for (uint i=0; i<triangle.size(); i++) {
					
					Geometry::Point &ca = triangle[i];
					Geometry::Point &cb = triangle[(i+1) % triangle.size()];
					cv::line(screen, ca, cb, cv::Scalar(64,64,64));
				}			
			}
		}

		if (Geometry::cut_polygon.size()>1) {
			
			for (size_t i = 0; i<Geometry::cut_polygon.size(); i++) {
				Geometry::Point &ca = Geometry::cut_polygon[i];
				Geometry::Point &cb = Geometry::cut_polygon[(i+1)%Geometry::cut_polygon.size()];
				cv::line(screen, ca, cb, cv::Scalar(64,64,255));


				//cv::Mat3b screen44;
				//cv::resize(screen, screen44, cv::Size(), 4, 4, cv::INTER_NEAREST);
				//cv::imshow("poligon", screen44);

				//cv::waitKey(1000);
			}
		}

		for (uint h=0; h<holes.size(); h++) {
			
			auto &hole = holes[h];
			
			if (hole.size()>1) {
				for (size_t i = 0; i<hole.size(); i++) {
					Geometry::Point &ca = hole[i];
					Geometry::Point &cb = hole[(i+1)%hole.size()];
					cv::line(screen, ca, cb, cv::Scalar(64,64,255));
				}
			}			
		}
		


		cv::Mat3b screen44;
		cv::resize(screen, screen44, cv::Size(), 4, 4, cv::INTER_NEAREST);
		cv::imshow("poligon", screen44);
	}
	
	
}
