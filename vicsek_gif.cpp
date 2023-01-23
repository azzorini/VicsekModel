#include <iostream>
#include <fstream>
#include <string>
#include <list>
#include <random>
#include <cmath>

const unsigned SEED = 2802;

class Bird {
	// Position
	double x, y;
	// Velocity
	double v;
	// Angle
	double theta;

	// Time step
	double h;

	// Lenght of the box in which the bird is
	double L;

	// Radius at which it can see other birds
	double r_sq;

	// Random numbers
	std::mt19937* gen;
	std::uniform_real_distribution<double> ran_theta;
	std::uniform_real_distribution<double> ran_u;

	// Variables to do a partial update
	double x_new, y_new, theta_new;
public:
	// Default constructor
	Bird(std::mt19937* g, double length = 7., double vel = .03, double eta = 2., double rad = 1, double dt = 1)
	: v(vel), h(dt), L(length), r_sq(rad*rad), gen(g), ran_theta(-eta/2, eta/2), ran_u(0, 1) {
		theta = 2*M_PI*ran_u(*gen);
		x = L*ran_u(*gen);
		y = L*ran_u(*gen);
	}

	// Copy constructor
	Bird(const Bird& o)
	: x(o.x), y(o.y), v(o.v), theta(o.theta), h(o.h), L(o.L), r_sq(o.r_sq), gen(o.gen),
	ran_theta(o.ran_theta), ran_u(o.ran_u), x_new(o.x_new), y_new(o.y_new), theta_new(o.theta_new) {}

	// Setters
	void set_pos(double x_, double y_ = -1) {
		if (x_ > -0.01) x = x_;
		if (y_ > -0.01) y = y_;
	}
	void set_angle(double angle) { theta = angle; }
	void set_v(double vel) { v = vel; }
	void set_dt(double dt) { h = dt; }
	void set_eta(double eta) { ran_theta = std::uniform_real_distribution<double>(-eta/2, eta/2); }
	void set_L(double, double, double);
	void set_r(double rad) { r_sq = rad*rad; }

	// Getters
	double get_x() const { return x; }
	double get_y() const { return y; }
	double get_v() const { return v; }
	double get_dt() const { return h; }
	double get_angle() const { return theta; }
	double get_L() const { return L; }
	double get_r() const { return std::sqrt(r_sq); }
	
	// Compute the square of the distance using periodic boundary conditions
	double distance_sq(const Bird&) const;

	// Computes the average angle fo the birds in the list
	// that are at a distance r from this bird
	double average_angle(const std::list<Bird>&) const;

	// Updates the angle and the position of the bird
	// it receives a list of all the birds to do so
	virtual void update(const std::list<Bird>&);

	// Same as update but instead of updating it
	// it saves everything in a new variable
	virtual void update_to_new(const std::list<Bird>&);

	// It loads the content of the new variables
	// into the usual variables
	void update_from_new() {
		x = x_new;
		y = y_new;
		theta = theta_new;
	}
};

// Sets the length of the box in which the bird is
// if the position exceeds the size of the new box
// the new position will be determined randomly
// for that reason it allows to give a position
void Bird::set_L(double length, double x_ = -1, double y_ = -1) {
	if (x > -0.01) x = x_;
	if (y > -0.01) y = y_;

	L = length;

	if (x > L) x = L*ran_u(*gen);
	if (y > L) y = L*ran_u(*gen);
}

double Bird::distance_sq(const Bird& b) const {
	double dx = std::abs(x - b.x);
	double dy = std::abs(y - b.y);

	if (dx > L/2) dx = L - dx;
	if (dy > L/2) dy = L - dy;

	return dx*dx + dy*dy;
}

double Bird::average_angle(const std::list<Bird>& birds) const {
	double sum_sin = .0, sum_cos = .0;

	for (auto it = birds.cbegin(); it != birds.cend(); it++) {
		if (distance_sq(*it) < r_sq) {
			sum_sin += std::sin(it->theta);
			sum_cos += std::cos(it->theta);
		}
	}

	// Technically what we should do is the average of the sums
	// but since std::atan2 only care about the proportion
	// there is no point in doing the averages
	return std::atan2(sum_sin, sum_cos);
}

void Bird::update(const std::list<Bird>& birds) {
	double avg_angle = average_angle(birds);

	// We update the position
	x += v*std::cos(theta)*h;
	y += v*std::sin(theta)*h;

	// We check that the bird has not crossed the borders
	// otherwise we put it in the right place
	if (x < 0) x += L;
	if (x > L) x -= L;

	if (y < 0) y += L;
	if (y > L) y -= L;

	// We update the angle
	theta = avg_angle + ran_theta(*gen);
}

void Bird::update_to_new(const std::list<Bird>& birds) {
	// We update the position
	x_new = x + v*std::cos(theta)*h;
	y_new = y + v*std::sin(theta)*h;
	
	// We check that the bird has not crossed the borders
	// otherwise we put it in the right place
	if (x_new < 0) x_new += L;
	if (x_new > L) x_new -= L;

	if (y_new < 0) y_new += L;
	if (y_new > L) y_new -= L;

	// We update the angle
	theta_new = average_angle(birds) + ran_theta(*gen);
}

int main() {
	// Simulation constants
	const unsigned N_birds = 300;
	const double vel = .03;
	const double L = 5.;
	const double eta = .1;
	const double r = 1.;
	const double dt = 1.;
	const unsigned N_steps = 5;
	const unsigned N = 400;

	std::mt19937 gen(SEED);
	std::list<Bird> bird_list;
	std::list<Bird>::iterator it;
	unsigned i, j;
	std::ofstream fbirds("birds_" + std::to_string(N_birds) + "_N_" + std::to_string(N)
			+ "_L_" + std::to_string(L) + "_movement_d.txt");

	for (i = 0; i < N_birds; i++) bird_list.push_back(Bird(&gen, L, vel, eta, r, dt));

	fbirds << "#t\tx\ty\ttheta\n";
	for (it = bird_list.begin(); it != bird_list.end(); it++)
		fbirds << 0 << '\t' << it->get_x() << '\t' << it->get_y() << '\t' << it->get_angle() << '\n';
	fbirds << '\n';
	for (i = 1; i < N; i++) {
		for (j = 0; j < N_steps; j++) {
			// We update every bird position to the new variables
			for (it = bird_list.begin(); it != bird_list.end(); it++) it->update_to_new(bird_list);

			// We load the new variables into the default variables and we save the data
			for (it = bird_list.begin(); it != bird_list.end(); it++) it->update_from_new();
		}
		for (it = bird_list.begin(); it != bird_list.end(); it++)
			fbirds << i*N_steps*dt << '\t' << it->get_x() << '\t' << it->get_y() << '\t' << it->get_angle() << '\n';
		fbirds << '\n';
	}

	fbirds.close();

	return 0;
}
