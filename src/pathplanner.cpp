#include <iostream>
#include <algorithm>
#include <math.h>
#include <time.h>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "pathplanner.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

#define M_PI       3.14159265358979323846   // pi

double deg2rad(double x)
{
	return x * M_PI / 180;
}

double rad2deg(double x)
{
	return x * 180 / M_PI;
}

double mph2ms(double mph)
{
	return mph * 1609.34 / 3600.00;
}

double ms2mph(double ms)
{
	return ms * 3600.00 / 1609.34;
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
}

void normalize_vector(double &nx, double &ny)
{
	double norm = sqrt(nx * nx + ny * ny);

	nx /= norm;
	ny /= norm;
}

vector<double> diff_poly(const vector<double> &coef)
{
	vector<double> diff_coef;

	int coef_num = coef.size();
	for (int i = 1; i < coef_num; ++i)
	{
		diff_coef.push_back(coef[i] * i);
	}

	return diff_coef;
}

double calc_poly(const vector<double> &coef, const double &T)
{
	double res = 0;

	double t = 1;
	int coef_num = coef.size();

	for (int i = 0; i < coef_num; ++i)
	{
		res += coef[i] * t;
		t *= T;
	}

	return res;
}

double calc_curvature(const vector<double> &d1, const vector<double> &d2, const double &T)
{
	return fabs(calc_poly(d2, T)) / pow(1.0 + pow(calc_poly(d1, T), 2.0), 1.5);
}

double logistics_function(const double &x)
{
	return 2.0 / (1 + exp(-x)) - 1.0;
}

// SegmentParameters

SegmentParameters::SegmentParameters()
{
}

SegmentParameters::SegmentParameters(const double &x, const double &y, const double &nx, const double &ny, const double &p, const double &s, const double &cv) : SegmentParameters()
{
	this->x = x;
	this->y = y;
	this->nx = nx;
	this->ny = ny;
	this->p = p;
	this->s = s;
	this->cv = cv;
}

SegmentParameters::SegmentParameters(const SegmentParameters &parameters) : SegmentParameters()
{
	copy_parameters(parameters);
}

SegmentParameters::~SegmentParameters()
{
}

SegmentParameters& SegmentParameters::operator=(const SegmentParameters &parameters)
{
	if (&parameters != this)
	{
		copy_parameters(parameters);
	}

	return *this;
}

void SegmentParameters::copy_parameters(const SegmentParameters &parameters)
{
	x = parameters.x;
	y = parameters.y;
	nx = parameters.nx;
	ny = parameters.ny;
	p = parameters.p;
	s = parameters.s;
	cv = parameters.cv;
}

// MapSegment

MapSegment::MapSegment()
{
}

MapSegment::MapSegment(const double &x0, const double &y0, const double &nx0, const double &ny0,
	const double &x1, const double &y1, const double &nx1, const double &ny1,
	const double &s0, const double &s1) : MapSegment()
{
	// Asume road map segment is non-null size and angle change less than 90 degree.

	this->s0 = s0;
	double ds = s1 - s0;

	double nx = nx0;
	double ny = ny0;
	normalize_vector(nx, ny);
	
	double p = calc_normal_eq_p(x0, y0, nx, ny);

	sp0 = new SegmentParameters(x0, y0, nx, ny, p, ds, 0);

	double lx1, ly1;
	from_global(x1, y1, lx1, ly1);

	double lnx1, lny1;
	from_global(x0 + nx1, y0 + ny1, lnx1, lny1);
	normalize_vector(lnx1, lny1);

	double lx1_2 = lx1 * lx1;
	double lx1_3 = lx1_2 * lx1;
	double dy1 = -lnx1 / lny1;

	spline = { 0, 0, (3 * ly1 - dy1 * lx1) / lx1_2, (dy1 * lx1 - 2 * ly1) / lx1_3 };
	dspline = diff_poly(spline);
	ddspline = diff_poly(dspline);

	spi.push_back(SegmentParameters(0, 0, 0, 1, 0, 0, calc_curvature(dspline, ddspline, 0)));
	double last_y = 0;
	double last_x = 0;
	double total_s = 0;
	for (double x = -delta; x > lx1; x -= delta)
	{
		double y = calc_poly(spline, x);
		double dy = calc_poly(dspline, x);
		double norm = sqrt(1 + dy * dy);
		double nx = -dy / norm;
		double ny = 1.0 / norm;
		double p = calc_normal_eq_p(x, y, nx, ny);

		total_s += sqrt((x - last_x) * (x - last_x) + (y - last_y) * (y - last_y));

		spi.push_back(SegmentParameters(x, y, nx, ny, p, total_s, calc_curvature(dspline, ddspline, (x + last_x) / 2.0)));

		last_x = x;
		last_y = y;
	}

	total_s += sqrt((lx1 - last_x) * (lx1 - last_x) + (ly1 - last_y) * (ly1 - last_y));

	spi.push_back(SegmentParameters(lx1, ly1, lnx1, lny1, calc_normal_eq_p(lx1, ly1, lnx1, lny1), total_s, calc_curvature(dspline, ddspline, (lx1 + last_x) / 2.0)));

	double norm_factor = ds / total_s;

	int segments_num = spi.size();
	for (int i = 0; i < segments_num; ++i)
	{
		spi[i].s *= norm_factor;
	}
}

MapSegment::MapSegment(const MapSegment &segment) : MapSegment()
{
	copy_segment(segment);
}

MapSegment::~MapSegment()
{
	if (sp0 != 0)
	{
		delete sp0;
		sp0 = 0;
	}
}

MapSegment& MapSegment::operator=(const MapSegment &segment)
{
	if (&segment != this)
	{
		copy_segment(segment);
	}
	return *this;
}

void MapSegment::copy_segment(const MapSegment &segment)
{
	s0 = segment.s0;

	if (sp0 != 0)
	{
		delete sp0;
		sp0 = 0;
	}
	
	sp0 = new SegmentParameters(*(segment.sp0));

	spline = segment.spline;
	dspline = segment.dspline;
	ddspline = segment.ddspline;
	spi = segment.spi;
}

double MapSegment::calc_normal_eq_p(const double &x, const double &y, const double &nx, const double &ny)
{
	return x * ny - y * nx;
}

void MapSegment::from_global(const double &x, const double &y, double &lx, double &ly)
{
	double dx = x - sp0->x;
	double dy = y - sp0->y;

	lx = dx * sp0->ny - dy * sp0->nx;
	ly = dx * sp0->nx + dy * sp0->ny;
}

void MapSegment::to_global(const double &lx, const double &ly, double &x, double &y)
{
	x = sp0->x + lx * sp0->ny + ly * sp0->nx;
	y = sp0->y + ly * sp0->ny - lx * sp0->nx;
}

double MapSegment::calc_line_distance(const SegmentParameters *parm, const double &x, const double &y)
{
	return calc_normal_eq_p(x, y, parm->nx, parm->ny) - parm->p;
}

void MapSegment::to_frenet(const double &x, const double &y, double &s, double &d)
{
	double lx, ly;
	from_global(x, y, lx, ly);
	
	int segments_num = spi.size();

	if (lx >= 0 || segments_num < 2)
	{
		s = s0;
		d = ly;

		return;
	}

	SegmentParameters* last_segment = &(spi[0]);
	double last_dist = calc_line_distance(last_segment, lx, ly);
	for (int i = 1; i < segments_num; ++i)
	{
		SegmentParameters* cur_segment = &(spi[i]);
		double cur_dist = calc_line_distance(cur_segment, lx, ly);

		if (last_dist <= 0 && cur_dist > 0)
		{
			d = (lx - last_segment->x) * last_segment->nx + (ly - last_segment->y) * last_segment->ny;

			double xd = last_segment->x + d * last_segment->nx;
			double yd = last_segment->y + d * last_segment->ny;
			double p1 = xd * last_segment->nx + yd * last_segment->ny;

			double norm = last_segment->nx * cur_segment->nx + last_segment->ny * cur_segment->ny;
			double xi = (p1 * cur_segment->nx + cur_segment->p * last_segment->ny) / norm;
			double yi = (p1 * cur_segment->ny - cur_segment->p * last_segment->nx) / norm;

			double xdi = xi - xd;
			double ydi = yi - yd;

			s = s0 + last_segment->s + (cur_segment->s - last_segment->s) * abs(last_dist) / sqrt(xdi * xdi + ydi * ydi);

			return;
		}

		last_segment = cur_segment;
		last_dist = cur_dist;
	}

	s = s0 + last_segment->s;
	d = (lx - last_segment->x) * last_segment->nx + (ly - last_segment->y) * last_segment->ny;
}

void MapSegment::from_frenet(const double &s, const double &d, double &x, double &y, double &cv)
{
	double lx, ly, lcv;
	int segments_num = spi.size();

	if (s <= s0 || segments_num < 2)
	{
		lx = 0;
		ly = d;
		lcv = 0;
	}
	else
	{
		double ls = s - s0;
		bool is_found = false;

		SegmentParameters* last_segment = &(spi[0]);
		for (int i = 1; i < segments_num; ++i)
		{
			SegmentParameters* cur_segment = &(spi[i]);

			if (ls >= last_segment->s && ls < cur_segment->s)
			{
				double xd = last_segment->x + d * last_segment->nx;
				double yd = last_segment->y + d * last_segment->ny;
				double p1 = xd * last_segment->nx + yd * last_segment->ny;

				double norm = last_segment->nx * cur_segment->nx + last_segment->ny * cur_segment->ny;
				double xi = (p1 * cur_segment->nx + cur_segment->p * last_segment->ny) / norm;
				double yi = (p1 * cur_segment->ny - cur_segment->p * last_segment->nx) / norm;

				double xdi = xi - xd;
				double ydi = yi - yd;

				double s_factor = (ls - last_segment->s) / (cur_segment->s - last_segment->s);

				lx = xd + s_factor * xdi;
				ly = yd + s_factor * ydi;
				lcv = last_segment->cv;

				is_found = true;

				break;
			}

			last_segment = cur_segment;
		}

		if (!is_found)
		{
			lx = last_segment->x + d * last_segment->nx;
			ly = last_segment->y + d * last_segment->ny;
			lcv = last_segment->cv;
		}
	}

	to_global(lx, ly, x, y);
	cv = lcv;
}

// PathPlanner

PathPlanner::PathPlanner()
{
}

PathPlanner::~PathPlanner()
{
}

vector<double> PathPlanner::frenet_vector(double s1, double d1, double s2, double d2)
{
	double ds = fmod(fmod(s2, map_max_s) + map_max_s, map_max_s) - fmod(fmod(s1, map_max_s) + map_max_s, map_max_s);
	double dd = d2 - d1;

	// Circular correction
	if (fabs(ds) > (map_max_s / 2.0))
	{
		ds += (ds > 0 ? -1.0 : 1.0) * map_max_s;
	}

	return { ds, dd };
}

double PathPlanner::frenet_distance(double s1, double d1, double s2, double d2)
{
	vector<double> vec = frenet_vector(s1, d1, s2, d2);
	return sqrt(vec[0] * vec[0] + vec[1] * vec[1]);
}

void PathPlanner::to_frenet(const double &x, const double &y, double &s, double &d)
{
	s = 0;
	d = 0;

	int segment_num = map_segments.size();
	int nearest_segment_idx = -1;
	double best_distance = 0;

	for (int i = 0; i < segment_num; ++i)
	{
		double cur_distance = distance(x, y, map_segments[i].sp0->x, map_segments[i].sp0->y);
		if (nearest_segment_idx < 0 || cur_distance < best_distance)
		{
			nearest_segment_idx = i;
			best_distance = cur_distance;
		}
	}

	if (nearest_segment_idx >= 0)
	{
		int cur_segment_idx = nearest_segment_idx;

		if (map_segments[nearest_segment_idx].calc_line_distance(map_segments[nearest_segment_idx].sp0, x, y) <= 0)
		{
			// Above segment normal, move forward

			int next_segment_idx = (cur_segment_idx + 1) % segment_num;
			while (map_segments[next_segment_idx].calc_line_distance(map_segments[next_segment_idx].sp0, x, y) <= 0)
			{
				cur_segment_idx = next_segment_idx;
				next_segment_idx = (cur_segment_idx + 1) % segment_num;

				if (cur_segment_idx == nearest_segment_idx)
				{
					break;
				}
			}
		}
		else
		{
			// Below segment normal, move backward

			int next_segment_idx = (cur_segment_idx - 1 + segment_num) % segment_num;
			while (map_segments[next_segment_idx].calc_line_distance(map_segments[next_segment_idx].sp0, x, y) > 0)
			{
				cur_segment_idx = next_segment_idx;
				next_segment_idx = (cur_segment_idx - 1 + segment_num) % segment_num;

				if (cur_segment_idx == nearest_segment_idx)
				{
					break;
				}
			}

			cur_segment_idx = next_segment_idx;
		}

		map_segments[cur_segment_idx].to_frenet(x, y, s, d);
	}
}

void PathPlanner::from_frenet(const double &s, const double &d, double &x, double &y)
{
	double cv;
	from_frenet(s, d, x, y, cv);
}

void PathPlanner::from_frenet(const double &s, const double &d, double &x, double &y, double &cv)
{
	x = 0;
	y = 0;
	cv = 0;

	double s_mod = fmod(fmod(s, map_max_s) + map_max_s, map_max_s);

	int segment_num = map_segments.size();
	for (int i = 0; i < segment_num; ++i)
	{
		if (map_segments[i].s0 <= s_mod && (map_segments[i].s0 + map_segments[i].sp0->s) > s_mod)
		{
			map_segments[i].from_frenet(s_mod, d, x, y, cv);

			break;
		}
	}
}

double PathPlanner::get_curvature(const double &s, const double &d)
{
	double x;
	double y;
	double cv;

	from_frenet(s, d, x, y, cv);

	return cv;
}

void PathPlanner::initialize_map(const vector<double> &map_x, const vector<double> &map_y, const vector<double> &map_s,
	const vector<double> &map_dx, const vector<double> &map_dy, const double &map_max_s)
{
	this->map_max_s = map_max_s;
	map_segments.clear();

	int map_size = map_x.size();

	if (map_size >= 2)
	{
		for (int i = 0; i < map_size - 1; ++i)
		{
			map_segments.push_back(MapSegment(map_x[i], map_y[i], map_dx[i], map_dy[i], map_x[i + 1], map_y[i + 1], map_dx[i + 1], map_dy[i + 1], map_s[i], map_s[i + 1]));
		}

		map_segments.push_back(MapSegment(map_x[map_size - 1], map_y[map_size - 1], map_dx[map_size - 1], map_dy[map_size - 1], map_x[0], map_y[0], map_dx[0], map_dy[0], map_s[map_size - 1], map_max_s));
	}
}

void PathPlanner::init(const vector<double> &map_x, const vector<double> &map_y, const vector<double> &map_s,
	const vector<double> &map_dx, const vector<double> &map_dy,
	const double &map_max_s, const double &lane_width, const int &lane_num, const double &refresh_time_sec)
{
	// Initialize path planner and validate input data.

	try
	{
		if (map_x.size() < 4 || map_y.size() < 4 || map_s.size() < 4 ||
			map_dx.size() < 4 || map_dy.size() < 4)
		{
			throw std::runtime_error("Map is empty or size is too small.");
		}

		if (map_max_s <= 0)
		{
			throw std::runtime_error("Map max length must be positive.");
		}

		this->lane_width = lane_width;
		if (this->lane_width <= 0)
		{
			throw std::runtime_error("Lane width must be positive.");
		}

		this->lane_num = lane_num;
		if (this->lane_num < 1)
		{
			throw std::runtime_error("Lane number must be at least 1.");
		}

		this->refresh_time_sec = refresh_time_sec;
		if (this->refresh_time_sec <= 0)
		{
			throw std::runtime_error("Refersh time must be positive.");
		}

		initialize_map(map_x, map_y, map_s, map_dx, map_dy, map_max_s);

		is_initialized = true;
	}
	catch (std::exception& ex)
	{
		is_initialized = false;
		cout << ex.what() << endl;
	}
}

// PathPlanner::Vehicle

PathPlanner::Vehicle::Vehicle()
{
}

PathPlanner::Vehicle::Vehicle(const PathPlanner* pathPlanner) : Vehicle()
{
	this->pathPlanner = (PathPlanner*)pathPlanner;
}

PathPlanner::Vehicle::Vehicle(const PathPlanner::Vehicle &vehicle) : Vehicle()
{
	copy_vehicle(vehicle);
}

PathPlanner::Vehicle& PathPlanner::Vehicle::operator=(const PathPlanner::Vehicle &vehicle)
{
	if (&vehicle != this)
	{
		copy_vehicle(vehicle);
	}

	return *this;
}

PathPlanner::Vehicle::~Vehicle()
{
}

void PathPlanner::Vehicle::copy_vehicle(const PathPlanner::Vehicle &vehicle)
{
	pathPlanner = vehicle.pathPlanner;
	id = vehicle.id;

	s = vehicle.s;
	s_v = vehicle.s_v;
	s_a = vehicle.s_a;
	d = vehicle.d;
	d_v = vehicle.d_v;
	d_a = vehicle.d_a;

	y_v = vehicle.y_v;
	y_a = vehicle.y_a;
}

void PathPlanner::Vehicle::process_forward(const double &t)
{
	s += s_v * t + 0.5 * s_a * t * t;
	d += d_v * t + 0.5 * d_a * t * t;

	double next_s_v = s_v + s_a * t;
	if (s_a < 0 &&
		next_s_v < speed_epsilon)
	{
		next_s_v = 0.5 * s_v;
	}

	s_v = next_s_v;
	d_v += d_a * t;
}

void PathPlanner::Vehicle::process_forward_cartesian(const double &x, const double &y, const double &yaw, const double &speed, const double &t, double &x1, double &y1) const
{
	x1 = x + speed * t * cos(yaw);
	y1 = y + speed * t * sin(yaw);
}

void PathPlanner::Vehicle::update_coordinates()
{
	// Calculate vehicle moving direction in frenet coordinates
	y_v = M_PI / 2.0;

	try
	{
		y_v = fmod(M_PI / 2.0 - atan2(d_v, s_v) + 2 * M_PI, 2 * M_PI);
	}
	catch (std::exception &ex)
	{
		y_v = M_PI / 2.0;
	}

	// Calculate vehicle acceleration direction in frenet coordinates
	y_a = M_PI / 2.0;

	try
	{
		y_a = fmod(M_PI / 2.0 - atan2(d_a, s_a) + 2 * M_PI, 2 * M_PI);
	}
	catch (std::exception &ex)
	{
		y_a = M_PI / 2.0;
	}
}

void PathPlanner::Vehicle::move_forward(const int &steps)
{
	if (steps >= 1)
	{
		for (int i = 0; i < steps; i++)
		{
			process_forward(pathPlanner->refresh_time_sec);
		}

		update_coordinates();
	}
}

void PathPlanner::Vehicle::move_forward_time(const double &time)
{
	double t = time;
	while (t > pathPlanner->refresh_time_sec)
	{
		process_forward(pathPlanner->refresh_time_sec);

		t -= pathPlanner->refresh_time_sec;
	}

	if (t > 0)
	{
		process_forward(t);
	}
}

void PathPlanner::Vehicle::to_cartesian(double &x, double &y, double &yaw, double &speed) const
{
	pathPlanner->from_frenet(s, d, x, y);

	Vehicle next_vehicle_state(*this);
	next_vehicle_state.process_forward(time_delta);

	double x1, y1;
	pathPlanner->from_frenet(next_vehicle_state.s, next_vehicle_state.d, x1, y1);

	double dx = x1 - x;
	double dy = y1 - y;

	yaw = 0;
	speed = 0;

	try
	{
		yaw = fmod(atan2(dy, dx) + 2 * M_PI, 2 * M_PI);
		speed = sqrt(dx * dx + dy * dy) / time_delta;
	}
	catch (std::exception &ex)
	{
		yaw = 0;
		speed = 0;
	}
}

void PathPlanner::Vehicle::from_cartesian(const double &x, const double &y, const double &yaw, const double &speed)
{
	pathPlanner->to_frenet(x, y, s, d);

	double x1, y1;
	process_forward_cartesian(x, y, yaw, speed, time_delta, x1, y1);

	double s1, d1;
	pathPlanner->to_frenet(x1, y1, s1, d1);

	vector<double> vec = pathPlanner->frenet_vector(s, d, s1, d1);
	s_v = vec[0] / time_delta;
	d_v = vec[1] / time_delta;

	// Initial acceleration cannot be calculated and can be ignored in this model
	s_a = 0;
	d_a = 0;
}

PathPlanner::Vehicle PathPlanner::create_vehicle_cartesian(const int &id, const double &x, const double &y, const double &yaw, const double &speed)
{
	Vehicle vehicle(this);

	vehicle.id = id;
	vehicle.from_cartesian(x, y, yaw, speed);
	vehicle.update_coordinates();

	return vehicle;
}

PathPlanner::Vehicle PathPlanner::create_vehicle_frenet(const int &id, const double &s, const double &s_v, const double &s_a, const double &d, const double &d_v, const double &d_a)
{
	Vehicle vehicle(this);

	vehicle.id = id;
	vehicle.s = s;
	vehicle.s_v = s_v;
	vehicle.s_a = s_a;
	vehicle.d = d;
	vehicle.d_v = d_v;
	vehicle.d_a = d_a;
	vehicle.update_coordinates();

	return vehicle;
}

// PathPlanner::PIDController
PathPlanner::PIDController::PIDController()
{
	reset();
}

PathPlanner::PIDController::PIDController(const PathPlanner* pathPlanner) : PIDController()
{
	this->pathPlanner = (PathPlanner*)pathPlanner;
}

PathPlanner::PIDController::PIDController(const PIDController &controller) : PIDController()
{
	copy_controller(controller);
}

PathPlanner::PIDController& PathPlanner::PIDController::operator=(const PIDController &controller)
{
	if (&controller != this)
	{
		copy_controller(controller);
	}

	return *this;
}

PathPlanner::PIDController::~PIDController()
{
}

void PathPlanner::PIDController::copy_controller(const PIDController &controller)
{
	pathPlanner = controller.pathPlanner;
	speed_s_history = controller.speed_s_history;
	distance_s_history = controller.distance_s_history;
	lateral_distance_points = controller.lateral_distance_points;
	speed_s_sum = controller.speed_s_sum;
	distance_s_sum = controller.distance_s_sum;
	lateral_distance_calc_speed = controller.lateral_distance_calc_speed;

	target_path_lane = controller.target_path_lane;
	target_speed = controller.target_speed;
	is_distance_limited = controller.is_distance_limited;
	target_distance = controller.target_distance;
}

void PathPlanner::PIDController::set_target_lane(const int &target_path_lane)
{
	if (this->target_path_lane != target_path_lane)
	{
		lateral_distance_points.clear();
	}

	this->target_path_lane = target_path_lane;
}

void PathPlanner::PIDController::set_target(const double &target_speed, const bool &is_distance_limited, const double &target_distance)
{
	this->target_speed = target_speed;
	this->is_distance_limited = is_distance_limited;
	this->target_distance = target_distance;

	if (!this->is_distance_limited)
	{
		distance_s_history.clear();
		distance_s_sum = 0;
	}
}

void PathPlanner::PIDController::reset()
{
	speed_s_history.clear();
	distance_s_history.clear();
	lateral_distance_points.clear();
	speed_s_sum = 0;
	distance_s_sum = 0;
	lateral_distance_calc_speed = 0;

	target_path_lane = -1;
	target_speed = 0;
	is_distance_limited = 0;
	target_distance = 0;
}

void PathPlanner::PIDController::process(Vehicle &vehicle)
{
	double curvature = pathPlanner->get_curvature(vehicle.s, vehicle.d);
	double curv_speed_factor = fabs(1.0 + curvature * vehicle.d);
	double acc_norm = vehicle.s_v * vehicle.s_v * curv_speed_factor * curv_speed_factor * curvature;

	Vehicle vehicle_proj(vehicle);
	vehicle_proj.move_forward(1);

	double speed_s_dv = (vehicle_proj.s_v - vehicle.s_v) * curv_speed_factor;
	double speed_s_da = vehicle_proj.s_a - vehicle.s_a;
	double speed = sqrt(vehicle.s_v * vehicle.s_v * curv_speed_factor * curv_speed_factor + vehicle.d_v * vehicle.d_v);
	double speed_diff = target_speed - speed;
	
	speed_s_sum += speed_diff;
	speed_s_history.push_back(speed_diff);
	if (speed_s_history.size() > speed_history_size)
	{
		speed_s_sum -= speed_s_history[0];
		speed_s_history.erase(speed_s_history.begin());
	}

	double acc_s_delta = Kps * speed_diff + Kiv * speed_s_sum - Kdsv * speed_s_dv - Kda * acc_norm;

	if (is_distance_limited && vehicle.s_v > critical_speed)
	{
		double s_dist_delta = pathPlanner->frenet_vector(vehicle.s, 0, target_distance, 0)[0];

		distance_s_sum += s_dist_delta;
		distance_s_history.push_back(s_dist_delta);
		if (distance_s_history.size() > distance_history_size)
		{
			distance_s_sum -= distance_s_history[0];
			distance_s_history.erase(distance_s_history.begin());
		}

		acc_s_delta += Kid * distance_s_sum;

		if (s_dist_delta < 0)
		{
			acc_s_delta += Kpts * s_dist_delta * vehicle.s_v;

			double next_s_dist_delta = pathPlanner->frenet_vector(vehicle_proj.s, 0, target_distance, 0)[0];
			double s_ddist_ddelta = s_dist_delta - next_s_dist_delta;
			acc_s_delta += Kdts * s_ddist_ddelta;
		}
	}

	vehicle.s_a += acc_s_delta / curv_speed_factor;

	if (vehicle.s_v >= min_speed_to_update_d)
	{
		double target_d_distance = pathPlanner->lane_width * (target_path_lane + 0.5);
		double d_dist_delta = target_d_distance - vehicle.d;

		bool is_updated = false;

		if (lateral_distance_points.size() > 0 && fabs(lateral_distance_calc_speed - vehicle.s_v) >= max_lateral_distance_calc_speed_diff)
		{
			lateral_distance_points.clear();
			lateral_distance_calc_speed = 0;
		}

		if (fabs(d_dist_delta) > max_lateral_dist_to_lane_center || lateral_distance_points.size() > 0)
		{
			if (lateral_distance_points.size() <= 0 &&
				vehicle.s_v > min_speed_to_update_d)
			{
				int estimated_periods_num = min(max_lateral_distance_points_at_once, calc_lane_change_periods(vehicle));
				vector<double> coef_d = pathPlanner->jerk_minimized_trajectory({ vehicle.d, vehicle.d_v, vehicle.d_a }, { target_d_distance, 0, 0 }, estimated_periods_num * pathPlanner->refresh_time_sec);

				vector<double> coef_dd = diff_poly(coef_d);
				vector<double> coef_ddd = diff_poly(coef_dd);

				double t = 0;
				vector<double> avg_acc;
				for (int i = 0; i < estimated_periods_num; ++i)
				{
					double next_d = calc_poly(coef_d, t);
					double next_d_v = calc_poly(coef_dd, t);
					double next_d_a = calc_poly(coef_ddd, t);

					lateral_distance_points.push_back({ next_d, next_d_v, next_d_a });

					t += pathPlanner->refresh_time_sec;
				}

				lateral_distance_calc_speed = vehicle.s_v;
			}

			if (lateral_distance_points.size() > 0)
			{
				vehicle.d = lateral_distance_points[0][0];
				vehicle.d_v = lateral_distance_points[0][1];
				vehicle.d_a = lateral_distance_points[0][2];

				lateral_distance_points.erase(lateral_distance_points.begin());

				is_updated = true;
			}
		}
		
		if (!is_updated)
		{
			double speed_d_dv = vehicle_proj.d_v - vehicle.d_v;
			double speed_d_da = vehicle_proj.d_a - vehicle.d_a;

			double next_d_dist_delta = target_d_distance - vehicle_proj.d;
			double d_ddist_ddelta = d_dist_delta - next_d_dist_delta;
			double acc_d_delta = Kpd * d_dist_delta - Kdtd * d_ddist_ddelta;

			vehicle.d_a += acc_d_delta;
		}
	}
}

bool PathPlanner::PIDController::can_change_lane(const Vehicle &vehicle)
{
	return
		target_path_lane < 0 ||
		fabs(vehicle.d - pathPlanner->lane_width * (target_path_lane + 0.5)) <= max_lateral_dist_to_lane_center;
}

int PathPlanner::PIDController::calc_lane_change_periods(const Vehicle &vehicle)
{
	double target_d_distance = pathPlanner->lane_width * (target_path_lane + 0.5);
	double d_dist_delta = target_d_distance - vehicle.d;
	double estimated_time = fabs(d_dist_delta) / fabs(vehicle.s_v * sin(deg2rad(change_lane_angle_deg)));
	int estimated_periods_num = estimated_time / pathPlanner->refresh_time_sec + 1.0;

	return estimated_periods_num;
}

// PathPlanner::PathStep

PathPlanner::PathStep::PathStep(const Vehicle &vehicle, const PIDController &pidController)
{
	this->vehicle.push_back(vehicle);
	this->pidController.push_back(pidController);
}

PathPlanner::PathStep::~PathStep()
{
}

const PathPlanner::Vehicle& PathPlanner::PathStep::getVehicle() const
{
	return vehicle[0];
}

const PathPlanner::PIDController& PathPlanner::PathStep::getPIDController() const
{
	return pidController[0];
}

// PathPlanner::RoadStatistic

PathPlanner::RoadStatistic::RoadStatistic()
{
}

PathPlanner::RoadStatistic::~RoadStatistic()
{
}

void PathPlanner::RoadStatistic::reset_lane_history()
{
	lane_avg_speed_history.clear();
	lane_sum_speed.clear();
}

void PathPlanner::RoadStatistic::push_lane_speed(const vector<double> &lane_speed)
{
	int lane_num = lane_speed.size();

	while (lane_sum_speed.size() < lane_num)
	{
		lane_sum_speed.push_back(0);
	}

	for (int i = 0; i < lane_num; ++i)
	{
		lane_sum_speed[i] += lane_speed[i];
	}

	lane_avg_speed_history.push_back(lane_speed);
	if (lane_avg_speed_history.size() > lane_avg_speed_history_size)
	{
		int min_lane_num = min(lane_num, int(lane_avg_speed_history[0].size()));

		for (int i = 0; i < min_lane_num; ++i)
		{
			lane_sum_speed[i] -= lane_avg_speed_history[0][i];
		}

		lane_avg_speed_history.erase(lane_avg_speed_history.begin());
	}
}

double PathPlanner::RoadStatistic::get_lane_avg_speed(const int &lane_num)
{
	if (lane_avg_speed_history.size() <= 0)
	{
		return 0;
	}

	if (lane_num >= lane_sum_speed.size())
	{
		return 0;
	}

	return lane_sum_speed[lane_num] / lane_avg_speed_history.size();
}

// PathPlanner

bool PathPlanner::check_collision(const vector<PathStep> &plan, const vector<Vehicle> &vehicles, const bool &check_safe_buffer)
{
	int plan_size = plan.size();
	int vehicles_num = vehicles.size();
	int check_full_collision_num = check_full_collision_sec / refresh_time_sec;

	if (plan_size > 0 && vehicles_num > 0)
	{
		double vehicle_min_size = min(vehicle_max_width, vehicle_max_length);

		vector<Vehicle> vehicles_next;
		for (int v = 0; v < vehicles_num; ++v)
		{
			vehicles_next.push_back(vehicles[v]);
		}

		vehicles_num = vehicles_next.size();

		for (int i = 0; i < plan_size; ++i)
		{
			const Vehicle &hero_pos = plan[i].getVehicle();
			double hero_d1 = hero_pos.d + vehicle_max_width / 2.0;
			double hero_d2 = hero_pos.d - vehicle_max_width / 2.0;

			if (i == check_full_collision_num)
			{
				vector<Vehicle> new_vehicles;
				for (int v = 0; v < vehicles_num; ++v)
				{
					// Vehicle is ahead
					if (frenet_vector(hero_pos.s, 0, vehicles_next[v].s, 0)[0] > 0)
					{
						new_vehicles.push_back(vehicles_next[v]);
					}
				}

				vehicles_next = new_vehicles;
				vehicles_num = vehicles_next.size();
				if (vehicles_num <= 0)
				{
					break;
				}
			}

			for (int v = 0; v < vehicles_num; ++v)
			{
				const Vehicle &vehicle_pos = vehicles_next[v];

				double vehicle_d1 = vehicle_pos.d + vehicle_max_width / 2.0;
				double vehicle_d2 = vehicle_pos.d - vehicle_max_width / 2.0;

				if (vehicle_d1 < hero_d2 && vehicle_d2 > hero_d1)
				{
					double vehicle_dist = fabs(frenet_vector(hero_pos.s, 0, vehicles_next[v].s, 0)[0]);
					if (check_safe_buffer &&
						vehicle_dist < (safe_buffer_collision + vehicle_max_length))
					{
						return false;
					}
				}

				// First check just distance between vehicle centers
				if (frenet_distance(hero_pos.s, hero_pos.d, vehicle_pos.s, vehicle_pos.d) <= vehicle_min_size)
				{
					return false;
				}

				// Next check distance by ellipse model
				double focus = sqrt(fabs(vehicle_max_length * vehicle_max_length - vehicle_max_width * vehicle_max_width));
				vector<double> f1 = { vehicle_pos.s + focus * sin(vehicle_pos.y_v), vehicle_pos.d + focus * cos(vehicle_pos.y_v) };
				vector<double> f2 = { vehicle_pos.s - focus * sin(vehicle_pos.y_v), vehicle_pos.d - focus * cos(vehicle_pos.y_v) };
				if (frenet_distance(hero_pos.s, hero_pos.d, f1[0], f1[1]) + frenet_distance(hero_pos.s, hero_pos.d, f2[0], f2[1]) < 2.0 * vehicle_max_length)
				{
					return false;
				}

				vehicles_next[v].move_forward(1);
			}
		}
	}

	return true;
}

vector<double> PathPlanner::jerk_minimized_trajectory(const vector< double> &start, const vector <double> &end, const double &T)
{
	vector<double> t;
	t.push_back(1);
	for (int i = 0; i < 5; ++i)
	{
	t.push_back(t[i] * T);
	}

	MatrixXd A(3, 3);
	A <<
	t[3], t[4], t[5],
	3 * t[2], 4 * t[3], 5 * t[4],
	6 * t[1], 12 * t[2], 20 * t[3];

	VectorXd y(3);
	y <<
	end[0] - start[0] - start[1] * t[1] - 0.5 * start[2] * t[2],
	end[1] - start[1] - start[2] * t[1],
	end[2] - start[2];

	VectorXd x = A.fullPivHouseholderQr().solve(y);

	return { start[0], start[1], 0.5 * start[2], x[0], x[1], x[2] };
}

vector<PathPlanner::PathStep> PathPlanner::generate_path(const PathStep &hero_car_status, const vector<Vehicle> &vehicles, const double &speed_limit)
{
	int vehicles_num = vehicles.size();
	
	// Vehicle ahead hero vehicle in appropriate lane
	vector<const Vehicle*> lane_vehicles_ahead;
	vector<double> lane_vehicles_ahead_distance;

	// Vehicle behind hero vehicle in appropriate lane
	vector<const Vehicle*> lane_vehicles_behind;
	vector<double> lane_vehicles_behind_distance;

	// If lane has vehicle near to hero vehicle which prevents lane change
	vector<bool> lane_collision;

	// Maximal distance used for planning
	const double speed_up_limit = speed_limit * keep_below_speed_limit_factor;
	const double horizon_distance = speed_limit * max_planned_time_sec;
	const double long_horizon_distance = long_horizon_distance_mult * horizon_distance;

	int hero_lane = -1;
	double hero_lane_center_distance = 0;

	const Vehicle &hero_car = hero_car_status.getVehicle();
	PIDController pid_upd(hero_car_status.getPIDController());

	// Minimal speed for vehicles ahead hero in the line
	// It will be used for average calculation
	vector<double> lane_min_speed;
	for (int i = 0; i < lane_num; ++i)
	{
		lane_min_speed.push_back(speed_up_limit);
	}

	double lane_start = 0;
	for (int i = 0; i < lane_num; ++i)
	{
		double lane_end = lane_start + lane_width;
		double lane_center = lane_start + lane_width / 2.0;
		bool has_collision = false;
			
		const Vehicle* vehicle_ahead = 0;
		double distance_to_ahead_vehicle = 0;

		const Vehicle* vehicle_behind = 0;
		double distance_to_behind_vehicle = 0;

		// Find vehicles just before, just after hero, and aside of hero
		// It block line and limit hero speed and maneuverability
		// It's possible that one vehicle blocks two lines simultaniously
		for (int v = 0; v < vehicles_num; ++v)
		{
			const Vehicle &vehicle_pos = vehicles[v];

			double d1 = vehicle_pos.d + vehicle_max_width / 2.0;
			double d2 = vehicle_pos.d - vehicle_max_width / 2.0;

			if (max(d1, d2) > lane_start && min(d1, d2) < lane_end)
			{
				bool include_in_avg_speed = false;
				double vehicle_distance = frenet_vector(hero_car.s, 0, vehicle_pos.s, 0)[0];

				if (vehicle_distance > vehicle_max_length)
				{
					if (vehicle_ahead == 0 || distance_to_ahead_vehicle > vehicle_distance)
					{
						distance_to_ahead_vehicle = vehicle_distance;
						vehicle_ahead = &(vehicles[v]);
					}

					include_in_avg_speed = true;
				}
				else if (vehicle_distance < vehicle_max_length)
				{
					if (vehicle_behind == 0 || distance_to_behind_vehicle > fabs(vehicle_distance))
					{
						distance_to_behind_vehicle = fabs(vehicle_distance);
						vehicle_behind = &(vehicles[v]);
					}
				}
				else
				{
					has_collision = true;
					include_in_avg_speed = true;
				}
				
				if (include_in_avg_speed && vehicle_distance <= long_horizon_distance)
				{
					lane_min_speed[i] = min(lane_min_speed[i], max(0.0, vehicle_pos.s_v));
				}
			}
		}

		if (vehicle_ahead != 0 && fabs(frenet_vector(hero_car.s, 0, vehicle_ahead->s, 0)[0]) > long_horizon_distance)
		{
			vehicle_ahead = 0;
		}

		if (vehicle_behind != 0 && fabs(frenet_vector(hero_car.s, 0, vehicle_behind->s, 0)[0]) > long_horizon_distance)
		{
			vehicle_behind = 0;
		}

		double cur_lane_center_distance = fabs(hero_car.d - lane_center);
		if (hero_lane < 0 || cur_lane_center_distance < hero_lane_center_distance)
		{
			hero_lane = i;
			hero_lane_center_distance = cur_lane_center_distance;
		}

		lane_vehicles_ahead.push_back((Vehicle*)vehicle_ahead);
		lane_vehicles_ahead_distance.push_back(distance_to_ahead_vehicle);

		lane_vehicles_behind.push_back((Vehicle*)vehicle_behind);
		lane_vehicles_behind_distance.push_back(distance_to_behind_vehicle);

		lane_collision.push_back(has_collision);

		lane_start = lane_end;
	}

	if (pid_upd.target_path_lane < 0)
	{
		pid_upd.set_target_lane(hero_lane);
	}

	bool pid_ready_change_lane = pid_upd.can_change_lane(hero_car);
	if (pid_ready_change_lane)
	{
		roadStatistic.push_lane_speed(lane_min_speed);
	}
	else
	{
		roadStatistic.reset_lane_history();
	}

	if (pid_ready_change_lane &&
		roadStatistic.lane_avg_speed_history.size() >= roadStatistic.lane_avg_speed_history_size &&
		(roadStatistic.get_lane_avg_speed(hero_lane) - hero_car.s_v) <= min_lane_speed_diff_to_change)
	{
		// Plan change line
		double best_lane_speed = 0;
		double best_lane_idx = -1;

		int dist = 2;
		while (true)
		{
			int lane_dist = dist / 2;
			int dist_sign = (dist % 2) == 0 ? -1 : 1;
			if (lane_dist > max_lanes_to_change)
			{
				break;
			}

			int next_lane = hero_lane + dist_sign * lane_dist;
			double hero_lane_speed_limit = roadStatistic.get_lane_avg_speed(hero_lane);
			double next_lane_speed_limit = roadStatistic.get_lane_avg_speed(next_lane);

			if (next_lane >= 0 && next_lane < lane_num &&
				fabs(hero_car.s_v - next_lane_speed_limit) <= max_speed_diff_to_lane_change &&
		  		((lane_vehicles_ahead[hero_lane] != 0 &&
				  frenet_vector(hero_car.s, 0, lane_vehicles_ahead[hero_lane]->s, 0)[0] <= long_horizon_distance &&
				  (lane_vehicles_ahead[next_lane] == 0 ||
				   (lane_vehicles_ahead[next_lane] != 0 &&
				    frenet_vector(hero_car.s, 0, lane_vehicles_ahead[next_lane]->s, 0)[0] > long_horizon_distance))) ||
				 (lane_vehicles_ahead[hero_lane] != 0 && lane_vehicles_ahead[next_lane] != 0 &&
				  (next_lane_speed_limit - hero_lane_speed_limit) >= min_speed_diff_to_lane_change)))
			{
				const Vehicle* vehicle_ahead_hero = lane_vehicles_ahead[hero_lane];

				bool can_change_line =
					vehicle_ahead_hero == 0 ||
					frenet_vector(hero_car.s, 0, vehicle_ahead_hero->s, 0)[0] >= (safe_buffer_lane_change + vehicle_max_length);

				if (can_change_line)
				{
					double lane_change_time = pid_upd.calc_lane_change_periods(hero_car) * refresh_time_sec;

					for (int i = 0; i < lane_num; ++i)
					{
						const Vehicle* vehicle_ahead = lane_vehicles_ahead[i];
						double vehicle_ahead_buffer_correction = 0;
						if (vehicle_ahead != 0)
						{
							vehicle_ahead_buffer_correction = max(0.0, hero_car.s_v * lane_change_time + 0.5 * hero_car.s_a * lane_change_time * lane_change_time - vehicle_ahead->s_v * lane_change_time);
						}

						const Vehicle* vehicle_behind = lane_vehicles_behind[i];
						double vehicle_behind_buffer_correction = 0;
						if (vehicle_behind != 0)
						{
							vehicle_behind_buffer_correction = max(0.0, vehicle_behind->s_v * lane_change_time - hero_car.s_v * lane_change_time - 0.5 * hero_car.s_a * lane_change_time * lane_change_time);
						}

						if (i >= min(next_lane, hero_lane) && i <= max(next_lane, hero_lane))
						{
							if (lane_collision[i])
							{
								can_change_line = false;
								break;
							}

							if (vehicle_ahead != 0)
							{
								if (frenet_vector(hero_car.s, 0, vehicle_ahead->s, 0)[0] < (safe_buffer_size + vehicle_ahead_buffer_correction + vehicle_max_length))
								{
									can_change_line = false;
									break;
								}
							}

							if (vehicle_behind != 0)
							{
								if (frenet_vector(vehicle_behind->s, 0, hero_car.s, 0)[0] < (safe_buffer_lane_change + vehicle_behind_buffer_correction + vehicle_max_length))
								{
									can_change_line = false;
									break;
								}
							}
						}
					}
				}

				if (can_change_line)
				{
					bool is_next_lane_no_ahead =
						lane_vehicles_ahead[next_lane] == 0 ||
						(lane_vehicles_ahead[next_lane] != 0 &&
						 frenet_vector(hero_car.s, 0, lane_vehicles_ahead[next_lane]->s, 0)[0] > long_horizon_distance);

					if (best_lane_idx < 0 ||
						best_lane_speed < next_lane_speed_limit ||
						is_next_lane_no_ahead)
					{
						best_lane_speed = next_lane_speed_limit;
						best_lane_idx = next_lane;

						if (is_next_lane_no_ahead)
						{
							break;
						}
					}
				}
			}

			++dist;
		}

		if (best_lane_idx >= 0)
		{
			// Change not more than 1 line at once
			if (best_lane_idx > hero_lane)
			{
				pid_upd.set_target_lane(hero_lane + 1);
			}
			else if (best_lane_idx < hero_lane)
			{
				pid_upd.set_target_lane(hero_lane - 1);
			}
		}
		else
		{
			pid_upd.set_target_lane(hero_lane);
		}
	}

	// Generate base plan
	Vehicle hero_car_upd(hero_car);
	
	vector<Vehicle> lane_vehicles_ahead_upd;
	for (int i = 0; i < lane_num; ++i)
	{
		if (lane_vehicles_ahead[i] != 0)
		{
			lane_vehicles_ahead_upd.push_back(Vehicle(*lane_vehicles_ahead[i]));
		}
		else
		{
			lane_vehicles_ahead_upd.push_back(Vehicle(this));
		}
	}

	vector<PathStep> new_plan;
	new_plan.push_back({ hero_car_upd, pid_upd });

	for (double t = refresh_time_sec; t <= max_planned_time_sec; t += refresh_time_sec)
	{
		bool is_s_distance_limited = false;
		double target_s_distance = hero_car_upd.s + horizon_distance;
		bool is_lane_changing = !pid_upd.can_change_lane(hero_car_upd);
		double cur_safe_buffer = is_lane_changing ? safe_buffer_size : safe_buffer_lane_change;
		double target_speed = speed_up_limit;

		double lane_start = 0;
		for (int i = 0; i < lane_num; ++i)
		{
			Vehicle* vehicle_ahead = (lane_vehicles_ahead[i] != 0) ? &lane_vehicles_ahead_upd[i] : 0;

			double lane_end = lane_start + lane_width;

			double d1 = hero_car_upd.d + vehicle_max_width / 2.0;
			double d2 = hero_car_upd.d - vehicle_max_width / 2.0;

			if (max(d1, d2) > lane_start && min(d1, d2) < lane_end)
			{
				if (vehicle_ahead != 0)
				{
					double cur_target_ahead_distance =
						hero_car_upd.s
						+ fabs(frenet_vector(hero_car_upd.s, 0, vehicle_ahead->s, 0)[0])
						- (cur_safe_buffer + vehicle_max_length);

					if (cur_target_ahead_distance < target_s_distance)
					{
						is_s_distance_limited = true;
						target_s_distance = cur_target_ahead_distance;
					}
				}

				if (vehicle_ahead != 0 &&
					fabs(frenet_vector(hero_car_upd.s, 0, vehicle_ahead->s, 0)[0]) < horizon_distance)
				{
					if (lane_min_speed[i] < target_speed)
					{
						target_speed = lane_min_speed[i];
					}
				}
			}

			if (vehicle_ahead != 0)
			{
				vehicle_ahead->move_forward(1);
			}

			lane_start = lane_end;
		}

		pid_upd.set_target(target_speed, is_s_distance_limited, target_s_distance);
		pid_upd.process(hero_car_upd);
		hero_car_upd.move_forward(1);

		new_plan.push_back({ hero_car_upd, pid_upd });
	}

	// If collision, try save current lane or switch to non-busy lane
	if (!this->check_collision(new_plan, vehicles, true))
	{
		vector<vector<PathStep>> strategy_list;

		const int strategy_keep_line = 1;
		const int strategy_left_line = 2;
		const int strategy_right_line = 3;
		const int strategy_keep_line_no_condition = 4;

		for (int strategy = strategy_keep_line; strategy <= strategy_keep_line_no_condition; ++strategy)
		{
			PIDController pid_upd(hero_car_status.getPIDController());
			pid_upd.reset();

			Vehicle hero_car_upd(hero_car);

			vector<Vehicle> lane_vehicles_ahead_upd;
			for (int i = 0; i < lane_num; ++i)
			{
				if (lane_vehicles_ahead[i] != 0)
				{
					lane_vehicles_ahead_upd.push_back(Vehicle(*lane_vehicles_ahead[i]));
				}
				else
				{
					lane_vehicles_ahead_upd.push_back(Vehicle(this));
				}
			}

			switch (strategy)
			{
			case strategy_right_line:
				if ((hero_lane + 1) >= lane_num || lane_collision[hero_lane + 1])
				{
					continue;
				}

				pid_upd.set_target_lane(hero_lane + 1);

				break;

			case strategy_left_line:
				if ((hero_lane - 1) < 0 || lane_collision[hero_lane - 1])
				{
					continue;
				}
				
				pid_upd.set_target_lane(hero_lane - 1);
				
				break;

			case strategy_keep_line:
				if (lane_collision[hero_lane])
				{
					continue;
				}

				pid_upd.set_target_lane(hero_lane);

				break;

			default:
				pid_upd.set_target_lane(hero_lane);
				break;
			}

			vector<PathStep> strategy_plan;
			strategy_plan.push_back({ hero_car_upd, pid_upd });

			for (double t = refresh_time_sec; t <= max_planned_time_sec; t += refresh_time_sec)
			{
				bool is_s_distance_limited = false;
				double target_s_distance = hero_car_upd.s + horizon_distance;
				double cur_safe_buffer = safe_buffer_size;
				double target_speed = speed_up_limit;

				double lane_start = 0;
				for (int i = 0; i < lane_num; ++i)
				{
					Vehicle* vehicle_ahead = (lane_vehicles_ahead[i] != 0) ? &lane_vehicles_ahead_upd[i] : 0;

					double lane_end = lane_start + lane_width;

					double d1 = hero_car_upd.d + vehicle_max_width / 2.0;
					double d2 = hero_car_upd.d - vehicle_max_width / 2.0;

					if (max(d1, d2) > lane_start && min(d1, d2) < lane_end)
					{
						if (vehicle_ahead != 0)
						{
							double cur_target_ahead_distance =
								hero_car_upd.s
								+ fabs(frenet_vector(hero_car_upd.s, 0, vehicle_ahead->s, 0)[0])
								- (cur_safe_buffer + vehicle_max_length);

							if (cur_target_ahead_distance < target_s_distance)
							{
								is_s_distance_limited = true;
								target_s_distance = cur_target_ahead_distance;
							}
						}

						if (vehicle_ahead != 0 &&
							fabs(frenet_vector(hero_car_upd.s, 0, vehicle_ahead->s, 0)[0]) < horizon_distance)
						{
							if (!lane_min_speed[i] < target_speed)
							{
								target_speed = lane_min_speed[i];
							}
						}
					}

					if (vehicle_ahead != 0)
					{
						vehicle_ahead->move_forward(1);
					}

					lane_start = lane_end;
				}

				pid_upd.set_target(target_speed, is_s_distance_limited, target_s_distance);
				pid_upd.process(hero_car_upd);
				hero_car_upd.move_forward(1);

				strategy_plan.push_back({ hero_car_upd, pid_upd });
			}

			if (strategy == strategy_keep_line ||
				strategy == strategy_keep_line_no_condition ||
				this->check_collision(strategy_plan, vehicles, true))
			{
				strategy_list.push_back(strategy_plan);
				break;
			}
		}

		if (strategy_list.size() > 0)
		{
			new_plan = strategy_list[strategy_list.size() - 1];
		}
	}

	return new_plan;
}

void PathPlanner::process(const double &car_x, const double &car_y, const double &car_s, const double &car_d, const double &car_yaw, const double &car_speed,
	const vector<double> &prev_path_x, const vector<double> &prev_path_y, const double &end_path_s, const double &end_path_d,
	const vector<vector<double>> &sensor_fusion, const double &speed_limit,
	vector<double> &next_x_vals, vector<double> &next_y_vals)
{
	// Validate if planner was initialized
	if (!is_initialized)
	{
		cout << "Path planner must be initialized first." << endl;
		return;
	}

	double timestamp = (double)clock() / CLOCKS_PER_SEC;
	double communication_lag = min_communication_daly_sec;

	if (last_timestamp <= 0)
	{
		last_timestamp = timestamp;
	}
	else
	{
		double dt = timestamp - last_timestamp;

		if (dt > max_communication_daly_sec)
		{
			communication_lag_history.clear();
			communication_lag_sum = 0;
		}
		else
		{
			communication_lag_history.push_back(dt);
			communication_lag_sum += dt;

			if (communication_lag_history.size() > communication_lag_history_size)
			{
				communication_lag_sum -= communication_lag_history[0];
				communication_lag_history.erase(communication_lag_history.begin());
			}

			communication_lag = max(min_communication_daly_sec, communication_lag_sum / communication_lag_history.size());
		}
	}

	last_timestamp = timestamp;

	//cout << "Lag: " << communication_lag << endl;

	// Convert inbound values into Vehicle class
	const Vehicle cur_hero_car = create_vehicle_cartesian(-1, car_x, car_y, fmod(deg2rad(car_yaw) + 2 * M_PI, 2 * M_PI), max(0.0, mph2ms(car_speed)));
	PathStep hero_car_status = PathStep(cur_hero_car, PIDController(this));
	
	vector<Vehicle> vehicles;
	int vehicles_num = sensor_fusion.size();
	for (int i = 0; i < vehicles_num; ++i)
	{
		double car_yaw = 0;
		double car_speed = 0;

		try
		{
			car_yaw = fmod(atan2(sensor_fusion[i][4], sensor_fusion[i][3]) + 2 * M_PI, 2 * M_PI);
			car_speed = sqrt(sensor_fusion[i][3] * sensor_fusion[i][3] + sensor_fusion[i][4] * sensor_fusion[i][4]);
		}
		catch (std::exception& ex)
		{
			car_yaw = 0;
			car_speed = 0;
		}

		Vehicle cur_car = create_vehicle_cartesian(int(sensor_fusion[i][0]), sensor_fusion[i][1], sensor_fusion[i][2], car_yaw, car_speed);
		vehicles.push_back(cur_car);
	}

	// Try match not-processed part of path returned by simulator and internal plan.
	// In case plan match current car position 
	int prev_path_size = min(prev_path_x.size(), prev_path_y.size());
	int plan_size = path_planned.size();
	bool cleanup_plan = true;

	if (plan_size > 0 && prev_path_size > 0)
	{
		double path_last_x, path_last_y, path_last_yaw, path_last_speed;
		path_planned[plan_size - 1].getVehicle().to_cartesian(path_last_x, path_last_y, path_last_yaw, path_last_speed);
		if (distance(prev_path_x[prev_path_size - 1], prev_path_y[prev_path_size - 1], path_last_x, path_last_y) < path_epsilon)
		{
			// Returned path match saved, find nearest point
			int nearest_idx = -1;
			double best_distance = 0;
			for (int i = 0; i < plan_size; ++i)
			{
				const Vehicle &cur_vehicle = path_planned[i].getVehicle();
				double cur_distance = frenet_distance(cur_vehicle.s, cur_vehicle.d, cur_hero_car.s, cur_hero_car.d);
				if (nearest_idx < 0 || cur_distance < best_distance)
				{
					nearest_idx = i;
					best_distance = cur_distance;
				}
			}

			if (nearest_idx >= 0)
			{
				bool is_first = true;
				double forward_dist = 0;

				do
				{
					if (is_first)
					{
						is_first = false;
					}
					else
					{
						if (nearest_idx > 0)
						{
							--nearest_idx;
						}
						else
						{
							forward_dist = 0;
							break;
						}
					}

					const Vehicle &cur_vehicle = path_planned[nearest_idx].getVehicle();
					vector<double> vec = frenet_vector(cur_vehicle.s, cur_vehicle.d, cur_hero_car.s, cur_hero_car.d);
					if (fabs(vec[1]) > max_plan_distance)
					{
						nearest_idx = -1;
						break;
					}

					forward_dist = vec[0] * (cur_hero_car.s_v >= 0 ? 1.0 : -1.0);
				} while (forward_dist < 0);

				if (nearest_idx >= 0 && nearest_idx < plan_size - 1)
				{
					// Calculate time shift
					const Vehicle &cur_vehicle = path_planned[nearest_idx].getVehicle();
					const Vehicle &cur_vehicle1 = path_planned[nearest_idx + 1].getVehicle();
					vector<double> vec = frenet_vector(cur_vehicle.s, 0, cur_vehicle1.s, 0);
					double time_shift = 0;

					try
					{
						if (fabs(vec[0]) > 0)
						{
							time_shift = refresh_time_sec * (1.0 - forward_dist / max(fabs(forward_dist), fabs(vec[0])));
						}
					}
					catch (std::exception& ex)
					{
						time_shift = 0;
					}

					// Move hero car status to the first plan step
					hero_car_status = path_planned[nearest_idx + 1];

					// Adjust other vehicles position
					if (time_shift > 0)
					{
						for (int v = 0; v < vehicles_num; ++v)
						{
							vehicles[v].move_forward_time(time_shift);
						}
					}

					// Save remain part of the plan
					path_planned.erase(path_planned.begin(), path_planned.begin() + nearest_idx + 1);
					plan_size = path_planned.size();

					int remain_points = min(plan_size, int(ceil(communication_lag / refresh_time_sec)));
					if (remain_points <= 0)
					{
						path_planned.clear();
					}
					else if (remain_points < plan_size)
					{
						path_planned.erase(path_planned.begin() + remain_points, path_planned.end());
					}

					// Check if current plan has collisions, if yes, cancel it
					// But remain several leading points, because it will be passed by vehicle due to communication lag
					/*if (!check_collision(path_planned, vehicles, true))
					{
						int remain_points = min(plan_size, int(ceil(communication_lag_sec / refresh_time_sec)));
						if (remain_points <= 0)
						{
							path_planned.clear();
						}
						else if (remain_points < plan_size)
						{
							path_planned.erase(path_planned.begin() + remain_points, path_planned.end());
						}
					}*/

					plan_size = path_planned.size();

					cleanup_plan = false;
				}
			}
		}
	}

	if (cleanup_plan)
	{
		path_planned.clear();
		plan_size = 0;
	}

	// If plan buffer was processed more than half or previous plan was cancelled, generate new plan steps
	if (plan_size <= 0.5 * max_planned_time_sec / refresh_time_sec)
	{
		// If previous plan remain, move vehicle state predictions to new state
		int move_forward_steps = 0;

		if (plan_size <= 0)
		{
			path_planned.push_back(hero_car_status);
			plan_size = 1;
		}
		else
		{
			move_forward_steps = plan_size - 1;
		}

		const double speed_limit_ms = mph2ms(speed_limit);

		// Several iterations may be required to fulfill prediction buffer
		while (plan_size < 0.75 * max_planned_time_sec / refresh_time_sec)
		{
			hero_car_status = path_planned[plan_size - 1];
			if (move_forward_steps > 0)
			{
				for (int v = 0; v < vehicles_num; ++v)
				{
					vehicles[v].move_forward(move_forward_steps);
				}

				move_forward_steps = 0;
			}

			vector<PathStep> best_path = generate_path(hero_car_status, vehicles, speed_limit_ms);
			int best_path_size = best_path.size();

			if (best_path_size > 1)
			{
				for (int i = 1; i < best_path_size; ++i)
				{
					path_planned.push_back(best_path[i]);
				}

				plan_size = path_planned.size();
				move_forward_steps = best_path_size - 1;
			}
			else
			{
				cout << "No path!" << endl;
				break;
			}

			break;
		}

		hero_car_status = path_planned[plan_size - 1];
	}

	double last_path_x, last_path_y;
	double last_path_s, last_path_d;

	for (int i = 0; i < plan_size; i++)
	{
		double path_x, path_y, path_yaw, path_speed;
		const Vehicle &cur_vehicle = path_planned[i].getVehicle();
		cur_vehicle.to_cartesian(path_x, path_y, path_yaw, path_speed);
		next_x_vals.push_back(path_x);
		next_y_vals.push_back(path_y);

		last_path_x = path_x;
		last_path_y = path_y;
		last_path_s = cur_vehicle.s;
		last_path_d = cur_vehicle.d;
	}
}
