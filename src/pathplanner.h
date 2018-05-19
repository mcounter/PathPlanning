#ifndef PATHPLANNER_H
#define PATHPLANNER_H

#include <vector>

using namespace std;

// Describe map segment parameters
class SegmentParameters
{
private:
	// Private constructor
	SegmentParameters();

	// Copy class instance
	void copy_parameters(const SegmentParameters &parameters);

public:
	// Point on the curve
	double x, y;

	// Normal vector and shift for normal linear equation
	double nx, ny, p;

	// Curve length
	double s;

	// Curvature (1/R)
	double cv;

	// Public constructor
	SegmentParameters(const double &x, const double &y, const double &nx, const double &ny, const double &p, const double &s, const double &cv);

	// Public copy constructor
	SegmentParameters(const SegmentParameters &parameters);

	// Destructor
	virtual ~SegmentParameters();

	// Copy operator
	SegmentParameters& operator=(const SegmentParameters &parameters);
};

// Describe map segment
class MapSegment
{
private:
	// Private constructor
	MapSegment();

	// Copy class instance
	void copy_segment(const MapSegment &segment);

	// Equation of the normal line based on the normal vector
	double calc_normal_eq_p(const double &x, const double &y, const double &nx, const double &ny);

public:
	// Step size for curve length integral equation
	// Bigger values cause non-smooth vehicle behaviour, smaller - decrease performance
	const double delta = 0.15;

	// Segment position in Frenet coordinates (assume d0 == 0)
	double s0;

	// Segment parameters in global coordinates
	SegmentParameters* sp0 = 0;

	// Curve coefficients (3rd order spline)
	vector<double> spline;

	// Curve differential coefficients
	vector<double> dspline;

	// Curve 2nd differential coefficients
	vector<double> ddspline;

	// Precalculated road sub-segment parameters, including initial and final vector in related coordinates.
	// For N segments, N+1 parameters in the vector.
	vector<SegmentParameters> spi;

	// Public constructor
	MapSegment(const double &x0, const double &y0, const double &nx0, const double &ny0,
		const double &x1, const double &y1, const double &nx1, const double &ny1,
		const double &s0, const double &s1);

	// Copy constructor
	MapSegment(const MapSegment &segment);

	// Desctructor
	virtual ~MapSegment();

	// Copy operator
	MapSegment& operator=(const MapSegment &segment);

	// Convert from global coordinates to segment local coordinates
	void from_global(const double &x, const double &y, double &lx, double &ly);
	
	// Convert from segment local coordinates to global coordinates
	void to_global(const double &lx, const double &ly, double &x, double &y);

	// Calc distance to the line
	double calc_line_distance(const SegmentParameters *parm, const double &x, const double &y);

	// Convert from Cartesian to Frenet coordinates
	void to_frenet(const double &x, const double &y, double &s, double &d);
	
	// Convert from Frenet to Cartesian coordinates
	void from_frenet(const double &s, const double &d, double &x, double &y, double &cv);
};

// Path planner main class
class PathPlanner
{
private:
	// Implements vehicle bicycle model in frenet coordinates
	class Vehicle
	{
	private:
		// Small time shift to convert speed vectors from different coordinates
		const double time_delta = 0.01;

		// Smallest speed value to prevent vehicle unpredictable behaviour near zero value, when PID controller suggests full stop
		const double speed_epsilon = 1.0;

		// Reference to parent class
		PathPlanner *pathPlanner;

		// Private constructor
		Vehicle();

		// Copy class instance
		void copy_vehicle(const Vehicle &vehicle);

		// Evaluate next vehicle position
		void process_forward(const double &t);

		// Evaluate next vehicle position in Cartesian coordinates
		void process_forward_cartesian(const double &x, const double &y, const double &yaw, const double &speed, const double &t, double &x1, double &y1) const;

	public:
		// Vehicle Id
		int id;

		// Vehicle longitudinal position on the road
		// Because track is cyclic, this value could be ambiguous
		double s;

		// Vehicle longitudinal speed
		double s_v;

		// Vehicle longitudinal acceleration
		double s_a;

		// Vehicle lateral position on the road
		double d;

		// Vehicle lateral speed
		double d_v;

		// Vehicle lateral acceleration
		double d_a;

		// Speed vector angle
		double y_v;

		// Acceleration vector angle
		double y_a;

		// Public constructor
		Vehicle(const PathPlanner *pathPlanner);

		// Copy constructor
		Vehicle(const Vehicle &vehicle);

		// Destructor
		virtual ~Vehicle();

		// Copy operator
		Vehicle& operator=(const Vehicle &vehicle);

		// When coordinates are updated, recalculate some dependant values
		void update_coordinates();

		// Move vehicle forward on some steps (20 ms each)
		void move_forward(const int &steps);

		// Move vehicle forward in time (split on 20 ms steps)
		void move_forward_time(const double &time);

		// Convert from Frenet to Cartesian coordinates
		void to_cartesian(double &x, double &y, double &yaw, double &speed) const;

		// Convert from Cartesian to Frenet coordinates
		void from_cartesian(const double &x, const double &y, const double &yaw, const double &speed);
	};

	// Implements mixed PID controller with path panning
	class PIDController
	{
	private:
		// PID controller parameters
		const double Kps = 0.0015;
		const double Kdsv = 0.6;
		const double Kda = 0.001;
		const double Kiv = 0.000003;
		const double Kid = 0.0000003;
		const double Kpts = 0.0018;
		const double Kdts = 0.003;

		const double Kpd = 0.4;
		const double Kdtd = 100.0;

		// Size of speed history buffer (for integral speed value calculation)
		const int speed_history_size = 150;
		
		// Size of distance history buffer (for integral distance to target value calculation)
		const int distance_history_size = 150;

		// Maximal lateral distance from lane center, when assumed vehicle keep the lane
		const double max_lateral_dist_to_lane_center = 0.25;

		// For lane change trajectory, maximal number of point calculated at once
		// Used for performance optimization to prevent calculation hang for a long time
		const int max_lateral_distance_points_at_once = 500;

		// Minimal speed when lateral distance can be updated (to prevent twitch on start)
		const double min_speed_to_update_d = 5.0;

		// Critical speed when vehicle further deceleration can cause unnecessary full stop
		const double critical_speed = 10.0;

		// Lateral distance is calculated for some longitudinal speed and when it change much must be recalculated.
		const double max_lateral_distance_calc_speed_diff = 5.0;

		// Average vehicle degree for lane change.
		// Used to estimate lane change time.
		const double change_lane_angle_deg = 4.0;

		// Speed history vector
		vector<double> speed_s_history;

		// Distance to target history vector
		vector<double> distance_s_history;

		// Pre-calculated lateral trajectory points
		vector<vector<double>> lateral_distance_points;

		// Longitudinal speed of vehicle at the moment of lateral trajectory was clculated
		double lateral_distance_calc_speed;

		// Integral speed difference
		double speed_s_sum = 0;

		// Integral distance to target difference
		double distance_s_sum = 0;

		// Reference to parent class
		PathPlanner *pathPlanner;

		// private constructor
		PIDController();

		// Copy class instance
		void copy_controller(const PIDController &controller);

	public:
		// Target planned path lane, -1 as initial condition
		int target_path_lane = -1;

		// Targtet vehicle speed
		double target_speed = 0;

		// Indicates if vehicle ahead exists close enougth
		bool is_distance_limited = false;

		// Distance current vehicle must cover to be on minimal safe distance from the vehicle ahead 
		double target_distance = 0;

		// Public constructor
		PIDController(const PathPlanner *pathPlanner);

		// Copy constructor
		PIDController(const PIDController &controller);

		// Destructor
		virtual ~PIDController();

		// Copy operator
		PIDController& operator=(const PIDController &controller);

		// Set target lane
		void set_target_lane(const int &target_path_lane);

		// Set target speed and distance
		void set_target(const double &target_speed, const bool &is_distance_limited, const double &target_distance);

		// Reset controller
		void reset();

		// PID controller entry point. Update vehicle acceleration before predict next vehicle path point
		void process(Vehicle &vehicle);

		// Indicates that vehicle keep target line and ready to next lane change movement
		bool can_change_lane(const Vehicle &vehicle);

		// Estimates lane change time in 20 ms periods.
		int calc_lane_change_periods(const Vehicle &vehicle);
	};

	// Recorder path step
	// Plan passed to the simulator is saved in memory to keep all parameetrs and smoothly plan next path points
	class PathStep
	{
	private:
		// Path point, its hero vehicle state (vector is used here to simplify code, it's 1 element length)
		vector<Vehicle> vehicle;

		// PID controller state (vector is used here to simplify code, it's 1 element length)
		vector<PIDController> pidController;
		
	public:
		// Public constructor
		PathStep(const Vehicle &vehicle, const PIDController &pidController);

		// Destructor
		~PathStep();

		// Return vehicle state
		const Vehicle& getVehicle() const;

		// Return PID controller state
		const PIDController& getPIDController() const;
	};

	// Save road statistic (speed history)
	class RoadStatistic
	{
	public:
		// Size of history buffer
		const double lane_avg_speed_history_size = 150;

		// History buffer
		vector<vector<double>> lane_avg_speed_history;

		// Integral speed
		vector<double> lane_sum_speed;

		// Public constructor
		RoadStatistic();

		// Destructor
		~RoadStatistic();

		// Cleanup lane history buffer
		void reset_lane_history();

		// Add current speed limits for each lane in the buffer
		void push_lane_speed(const vector<double> &lane_speed);

		// Calculate average lane speed limit
		double get_lane_avg_speed(const int &lane_num);
	};

	// Maximum planned time
	const double max_planned_time_sec = 2.0;

	// Keep below speed limit factor
	const double keep_below_speed_limit_factor = 0.925;

	// Maximum acceptable distance to original plan
	const double max_plan_distance = 1.0;

	// Maximum vehicle width
	const double vehicle_max_width = 2.0;

	// Maximum vehicle length
	const double vehicle_max_length = 4.0;

	// Safe buffer to vehicle ahead for speed up to 50 MPH
	const double safe_buffer_size = 20.0;

	// Min safe buffer to vehicles safe to change line
	// This is a part of trajectory control for each point in each moment of time
	const double safe_buffer_lane_change = safe_buffer_size * 0.5;

	// Min safe buffer to recalculate path
	const double safe_buffer_collision = 5.0;

	// Minimal speed difference below current lane speed, when possible change the lane.
	const double min_lane_speed_diff_to_change = 2.5;

	// Minimum allowed error to match path points
	const double path_epsilon = 0.002;

	// Check full collision sec
	const double check_full_collision_sec = 0.2;

	// Maximal number of lanes to be changed at once
	// It's just for planning, only 1 lane can be changed at once
	const int max_lanes_to_change = 2;

	// Minimal speed difference to change lane, m/s
	// Otherwise it's assumed no sense to do this movement
	const double min_speed_diff_to_lane_change = 3.0;

	// Maximal speed difference to change line, m/s
	// It's assumed vehicle must speedup (or slowdown) first to change lane safly
	const double max_speed_diff_to_lane_change = 15.0;

	// Multiplier of when distance to vehicle ahead is ignored
	const double long_horizon_distance_mult = 1.5;

	// Communication delay
	const double min_communication_daly_sec = 0.02;

	// Maximal communication delay reset this stattistic
	const double max_communication_daly_sec = 5.0;

	// Size of the communication lag buffer
	const int communication_lag_history_size = 50;

	// Indicates that class is initialized successfully
	bool is_initialized = false;

	// Length of the track
	double map_max_s;

	// List of pre-processed map segments
	vector<MapSegment> map_segments;
	
	// Width of the one road lane.
	double lane_width;

	// Number of road lanes on hero vehicle side
	int lane_num;

	// Time distance between path adjacent points.
	double refresh_time_sec;

	// Last planned path
	vector<PathStep> path_planned;

	// Road statistic class instance
	RoadStatistic roadStatistic;

	// Last request time stamp
	double last_timestamp = 0;

	// History of total communication and calculation lags (time distance between requests)
	vector<double> communication_lag_history;

	// Integral communication lag
	double communication_lag_sum = 0;

	// Creates Frenet vector by 2 points
	// Right way to calculate longitudinal distance between 2 points in circular Frenet coordinates
	vector<double> frenet_vector(double s1, double d1, double s2, double d2);

	// Distance between 2 points in Frenet coordinates
	double frenet_distance(double s1, double d1, double s2, double d2);

	// Create map segments from loaded map points
	void initialize_map(const vector<double> &map_x, const vector<double> &map_y, const vector<double> &map_s,
		const vector<double> &map_dx, const vector<double> &map_dy, const double &map_max_s);

	// Convert Cartesian coorditanes to Frenet coordinates
	void to_frenet(const double &x, const double &y, double &s, double &d);
	
	// Convert Frenet coordinates to Cartesian coorditanes
	void from_frenet(const double &s, const double &d, double &x, double &y);

	// Convert Frenet coordinates to Cartesian coorditanes + return curvature
	void from_frenet(const double &s, const double &d, double &x, double &y, double &cv);

	// Return curvature of the road in some frenet point
	double get_curvature(const double &s, const double &d);

	// Create Vehicle class instance from data in Cartesian coordinates
	Vehicle create_vehicle_cartesian(const int &id, const double &x, const double &y, const double &yaw, const double &speed);

	// Create Vehicle class instance from data in Frenet coordinates
	Vehicle create_vehicle_frenet(const int &id, const double &s, const double &s_v, const double &s_a, const double &d, const double &d_v, const double &d_a);

	// Check planned path for collisisons
	bool check_collision(const vector<PathStep> &plan, const vector<Vehicle> &vehicles, const bool &check_safe_buffer = false);

	// Generates jerk minimized trajectory based on quintic polynomial
	vector<double> jerk_minimized_trajectory(const vector< double> &start, const vector <double> &end, const double &T);

	// Main path planner method
	vector<PathStep> generate_path(const PathStep &hero_car_status, const vector<Vehicle> &vehicles, const double &speed_limit);

public:
	
	// Public constructor
	PathPlanner();

	// Destructor
	virtual ~PathPlanner();

	// Initialize path planner instance with parameters and map
	virtual void init(const vector<double> &map_x, const vector<double> &map_y, const vector<double> &map_s,
		const vector<double> &map_dx, const vector<double> &map_dy,
		const double &map_max_s, const double &lane_width, const int &lane_num, const double &refresh_time_sec);

	// Main path planner entry point
	// Called from web server onMessage event handler
	virtual void process(const double &car_x, const double &car_y, const double &car_s, const double &car_d, const double &car_yaw, const double &car_speed,
		const vector<double> &prev_path_x, const vector<double> &prev_path_y, const double &end_path_s, const double &end_path_d,
		const vector<vector<double>> &sensor_fusion, const double &speed_limit,
		vector<double> &next_x_vals, vector<double> &next_y_vals);
};

#endif /*PATHPLANNER_H*/