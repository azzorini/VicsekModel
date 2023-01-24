#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <list>
#include <algorithm>
#include <random>
#include <cmath>

const unsigned SEED = 2802;

// Get the position in the 1D vector given by the indeces
inline unsigned pos_by_indeces(unsigned i, unsigned j, unsigned N_cells) {
	return i*N_cells + j;
}

class Bird {
	// Position
	double x, y;
	// Velocity
	double v;
	// Angle
	double theta;
	
	// Lenght of the box in which the bird is
	double L;

	// Link cell indeces
	unsigned i, j;
	// Variable that tell us if the numbers of cells exactly match the length
	bool exact_cells;

	// Random numbers
	std::mt19937* gen;
	std::uniform_real_distribution<double> ran_theta;
	std::uniform_real_distribution<double> ran_u;

	// Variables to do a partial update
	double x_new, y_new, theta_new;
public:
	// Default constructor
	Bird(std::mt19937*, double, double, double);

	// Copy constructor
	Bird(const Bird& o)
	: x(o.x), y(o.y), v(o.v), theta(o.theta), L(o.L), i(o.i), j(o.j), exact_cells(o.exact_cells), gen(o.gen),
	ran_theta(o.ran_theta), ran_u(o.ran_u), x_new(o.x_new), y_new(o.y_new), theta_new(o.theta_new) {}

	// Setters
	void set_pos(double x_, double y_ = -1) {
		if (x_ > -0.01) x = x_;
		if (y_ > -0.01) y = y_;
	}
	void set_angle(double angle) { theta = angle; }
	void set_v(double vel) { v = vel; }
	void set_eta(double eta) { ran_theta = std::uniform_real_distribution<double>(-eta/2, eta/2); }
	void set_L(double, double, double);

	// Getters
	double get_x() const { return x; }
	double get_y() const { return y; }
	double get_v() const { return v; }
	double get_angle() const { return theta; }
	double get_L() const { return L; }
	unsigned get_i() const { return i; }
	unsigned get_j() const { return j; }
	
	// Compute the square of the distance using periodic boundary conditions
	double distance_sq(const Bird&) const;

	// Computes the average angle fo the birds in the list
	// that are at a distance r from this bird
	double average_angle(const std::list<Bird>&, const std::vector< std::vector<unsigned> >&,
			const std::vector< std::list<Bird*> >&, unsigned) const;

	// Updates the angle and the position of the bird
	// it receives a list of all the birds to do so
	virtual void update(const std::list<Bird>&, const std::vector< std::vector<unsigned> >&,
			std::vector< std::list<Bird*> >&, unsigned);

	// Same as update but instead of updating it
	// it saves everything in a new variable
	virtual void update_to_new(const std::list<Bird>&, const std::vector< std::vector<unsigned> >&,
			const std::vector< std::list<Bird*> >&, unsigned);

	// It loads the content of the new variables
	// into the usual variables
	void update_from_new(std::vector< std::list<Bird*> >&, unsigned);

};

Bird::Bird(std::mt19937* g, double length = 7., double vel = .03, double eta = 2.)
: v(vel), L(length), exact_cells(!(std::ceil(length) > length)), gen(g), ran_theta(-eta/2, eta/2), ran_u(0, 1), x_new(0), y_new(0), theta_new(0) {
	theta = 2*M_PI*ran_u(*gen);
	x = L*ran_u(*gen);
	y = L*ran_u(*gen);
	i = (unsigned) x;
	j = (unsigned) y;
}

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

	if (dx > .5*L) dx = L - dx;
	if (dy > .5*L) dy = L - dy;

	return dx*dx + dy*dy;
}

double Bird::average_angle(const std::list<Bird>& birds, const std::vector< std::vector<unsigned> >& neighbours,
			const std::vector< std::list<Bird*> >& link_list, unsigned ind) const {
	double sum_sin = .0, sum_cos = .0;
	unsigned n, N_cells, N_neigh = neighbours.size();

	for (Bird* p : link_list[ind]) {
		if (distance_sq(*p) < 1) {
			sum_sin += std::sin(p->theta);
			sum_cos += std::cos(p->theta);
		}
	}

	for (n = 0; n < N_neigh; n++) {
		for (Bird* p : link_list[neighbours[n][ind]]) {
			if (distance_sq(*p) < 1) {
				sum_sin += std::sin(p->theta);
				sum_cos += std::cos(p->theta);
			}
		}
	}
	
	// Some additional checks if the number of cells does not exactly match the length
	if (!exact_cells) {
		N_cells = std::ceil(L);
		if (i == 0) {
			n = pos_by_indeces(N_cells-1, j, N_cells);
			for (Bird* p : link_list[neighbours[1][n]]) {
				if (distance_sq(*p) < 1) {
					sum_sin += std::sin(p->theta);
					sum_cos += std::cos(p->theta);
				}
			}
			for (Bird* p : link_list[neighbours[4][n]]) {
				if (distance_sq(*p) < 1) {
					sum_sin += std::sin(p->theta);
					sum_cos += std::cos(p->theta);
				}
			}
			for (Bird* p : link_list[neighbours[5][n]]) {
				if (distance_sq(*p) < 1) {
					sum_sin += std::sin(p->theta);
					sum_cos += std::cos(p->theta);
				}
			}
		} else if (i == N_cells - 2) {
				n = pos_by_indeces(N_cells-1, j, N_cells);
				for (Bird* p : link_list[neighbours[3][n]]) {
					if (distance_sq(*p) < 1) {
						sum_sin += std::sin(p->theta);
						sum_cos += std::cos(p->theta);
					}
				}
				for (Bird* p : link_list[neighbours[6][n]]) {
					if (distance_sq(*p) < 1) {
						sum_sin += std::sin(p->theta);
						sum_cos += std::cos(p->theta);
					}
				}
				for (Bird* p : link_list[neighbours[7][n]]) {
					if (distance_sq(*p) < 1) {
						sum_sin += std::sin(p->theta);
						sum_cos += std::cos(p->theta);
					}
				}

		}
		if (j == 0) {
			n = pos_by_indeces(i, N_cells-1, N_cells);
			for (Bird* p : link_list[neighbours[0][n]]) {
				if (distance_sq(*p) < 1) {
					sum_sin += std::sin(p->theta);
					sum_cos += std::cos(p->theta);
				}
			}
			for (Bird* p : link_list[neighbours[4][n]]) {
				if (distance_sq(*p) < 1) {
					sum_sin += std::sin(p->theta);
					sum_cos += std::cos(p->theta);
				}
			}
			for (Bird* p : link_list[neighbours[7][n]]) {
				if (distance_sq(*p) < 1) {
					sum_sin += std::sin(p->theta);
					sum_cos += std::cos(p->theta);
				}
			}
		} else if (j == N_cells - 2) {
				n = pos_by_indeces(i, N_cells-1, N_cells);
				for (Bird* p : link_list[neighbours[2][n]]) {
					if (distance_sq(*p) < 1) {
						sum_sin += std::sin(p->theta);
						sum_cos += std::cos(p->theta);
					}
				}
				for (Bird* p : link_list[neighbours[5][n]]) {
					if (distance_sq(*p) < 1) {
						sum_sin += std::sin(p->theta);
						sum_cos += std::cos(p->theta);
					}
				}
				for (Bird* p : link_list[neighbours[6][n]]) {
					if (distance_sq(*p) < 1) {
						sum_sin += std::sin(p->theta);
						sum_cos += std::cos(p->theta);
					}
				}

		}

	}

	// Technically what we should do is the average of the sums
	// but since std::atan2 only care about the proportion
	// there is no point in doing the averages
	return std::atan2(sum_sin, sum_cos);
}

void Bird::update(const std::list<Bird>& birds, const std::vector< std::vector<unsigned> >& neighbours,
			std::vector< std::list<Bird*> >& link_list, unsigned N_cells) {
	unsigned ind = pos_by_indeces(i, j, N_cells);
	bool changed = false;
	double avg_angle = average_angle(birds, neighbours, link_list, ind);
	unsigned aux;
	std::list<Bird*>::iterator it;

	// We update the position
	x += v*std::cos(theta);
	y += v*std::sin(theta);

	// We check that the bird has not crossed the borders
	// otherwise we put it in the right place
	if (x < 0) x += L;
	if (x > L) x -= L;

	if (y < 0) y += L;
	if (y > L) y -= L;

	// If we need to we update i and j
	aux = (unsigned) x;
	if (i != aux) {
		i = aux;
		changed = true;
	}

	aux = (unsigned) y;
	if (j != aux) {
		j = aux;
		changed = true;
	}

	if (changed) {
		// We update the link_list

		// We assume that this list is ok and we did find it
		it = std::find(link_list[ind].begin(), link_list[ind].end(), this);
		link_list[pos_by_indeces(i, j, N_cells)].push_back(*it);
		link_list[ind].erase(it);

	}

	// We update the angle
	theta = avg_angle + ran_theta(*gen);
}

void Bird::update_to_new(const std::list<Bird>& birds, const std::vector< std::vector<unsigned> >& neighbours,
			const std::vector< std::list<Bird*> >& link_list, unsigned N_cells) {
	// We update the position
	x_new = x + v*std::cos(theta);
	y_new = y + v*std::sin(theta);
	
	// We check that the bird has not crossed the borders
	// otherwise we put it in the right place
	if (x_new < 0) x_new += L;
	if (x_new > L) x_new -= L;

	if (y_new < 0) y_new += L;
	if (y_new > L) y_new -= L;

	// We update the angle
	theta_new = average_angle(birds, neighbours, link_list, pos_by_indeces(i, j, N_cells)) + ran_theta(*gen);
}

void Bird::update_from_new(std::vector< std::list<Bird*> >& link_list, unsigned N_cells) {
	unsigned aux, ind = pos_by_indeces(i, j, N_cells);
	bool changed = false;
	std::list<Bird*>::iterator it;

	x = x_new;
	y = y_new;
	theta = theta_new;

	// If we need to we update i and j
	aux = (unsigned) x;
	if (i != aux) {
		i = aux;
		changed = true;
	}

	aux = (unsigned) y;
	if (j != aux) {
		j = aux;
		changed = true;
	}

	if (changed) {
		// We update the link_list

		// We assume that this list is ok and we did find it
		it = std::find(link_list[ind].begin(), link_list[ind].end(), this);
		link_list[pos_by_indeces(i, j, N_cells)].push_back(*it);
		link_list[ind].erase(it);
	}
}

void compute_neighbours(std::vector< std::vector<unsigned> >& n, unsigned N_cells) {
	unsigned i, j;
	unsigned ant_i, next_i, ant_j, next_j;

	for (i = 0; i < N_cells; i++) {
		ant_i = (i > 0 ? i-1 : N_cells-1);
		next_i = (i < N_cells-1 ? i+1 : 0);
		for (j = 0; j < N_cells; j++) {
			ant_j = (j > 0 ? j-1 : N_cells-1);
			next_j = (j < N_cells-1 ? j+1 : 0);
			n[0][pos_by_indeces(i, j, N_cells)] = pos_by_indeces(ant_i, j, N_cells);
			n[1][pos_by_indeces(i, j, N_cells)] = pos_by_indeces(i, ant_j, N_cells);
			n[2][pos_by_indeces(i, j, N_cells)] = pos_by_indeces(next_i, j, N_cells);
			n[3][pos_by_indeces(i, j, N_cells)] = pos_by_indeces(i, next_j, N_cells);
			n[4][pos_by_indeces(i, j, N_cells)] = pos_by_indeces(ant_i, ant_j, N_cells);
			n[5][pos_by_indeces(i, j, N_cells)] = pos_by_indeces(next_i, ant_j, N_cells);
			n[6][pos_by_indeces(i, j, N_cells)] = pos_by_indeces(next_i, next_j, N_cells);
			n[7][pos_by_indeces(i, j, N_cells)] = pos_by_indeces(ant_i, next_j, N_cells);
		}
	}
}

void compute_link_list(std::list<Bird>& L, std::vector< std::list<Bird*> >& v, unsigned N_cells) {
	for (auto it = L.begin(); it != L.end(); it++)
		v[pos_by_indeces(it->get_i(), it->get_j(), N_cells)].push_back(&(*it));
}

double compute_va(const std::list<Bird>& L, unsigned N_birds) {
	double sum_sin = .0, sum_cos = .0;
	double ang;

	for (auto it = L.cbegin(); it != L.cend(); it++) {
		ang = it->get_angle();
		sum_sin += std::sin(ang);
		sum_cos += std::cos(ang);
	}

	return std::sqrt(sum_sin*sum_sin + sum_cos*sum_cos)/N_birds;
}

int main() {
	// Simulation constants
	const unsigned N_birds = 10000;
	const double rho = 4.;
	const double vel = .03;
	const double L = std::sqrt(N_birds/rho);
	const unsigned N_cells = (unsigned) std::ceil(L);
	const unsigned N_relax_ini = 1000;
	const unsigned N_relax = 500;
	const unsigned N_steps = 5;
	const unsigned N_measures = 1000;

	// Random number generator
	std::mt19937 gen(SEED);

	std::list<Bird> bird_list;
	std::list<Bird>::iterator it;
	std::vector< std::list<Bird*> > link_list(N_cells*N_cells, std::list<Bird*>());
	std::vector< std::vector<unsigned> > neighbours(8, std::vector<unsigned>(N_cells*N_cells));

	unsigned i, j;
	double va, va_avg, va_sq_avg;
	double eta = 5.;
	std::ofstream fout("VAvsETA_N_" + std::to_string(N_birds) + "_L_" + std::to_string(L) + "_i.txt");
	fout << "#eta\tv_a\tv_a err\n";

	compute_neighbours(neighbours, N_cells);

	for (i = 0; i < N_birds; i++) bird_list.push_back(Bird(&gen, L, vel, eta));

	compute_link_list(bird_list, link_list, N_cells);

	for (i = 0; i < N_relax_ini; i++) {
			// We update every bird position to the new variables
			for (it = bird_list.begin(); it != bird_list.end(); it++) it->update_to_new(bird_list, neighbours, link_list, N_cells);

			// We load the new variables into the default variables and we save the data
			for (it = bird_list.begin(); it != bird_list.end(); it++) it->update_from_new(link_list, N_cells);
	}

	while (eta > .05) {
		std::cout << "Starting calculation for eta = " << eta << "...\n";
		va_avg = va_sq_avg = .0;
		for (i = 0; i < N_relax; i++) {
				// We update every bird position to the new variables
				for (it = bird_list.begin(); it != bird_list.end(); it++) it->update_to_new(bird_list, neighbours, link_list, N_cells);

				// We load the new variables into the default variables and we save the data
				for (it = bird_list.begin(); it != bird_list.end(); it++) it->update_from_new(link_list, N_cells);
		}

		for (i = 0; i < N_measures; i++) {
			for (j = 0; j < N_steps; j++) {
				// We update every bird position to the new variables
				for (it = bird_list.begin(); it != bird_list.end(); it++) it->update_to_new(bird_list, neighbours, link_list, N_cells);

				// We load the new variables into the default variables and we save the data
				for (it = bird_list.begin(); it != bird_list.end(); it++) it->update_from_new(link_list, N_cells);
			}

			// We do the measure
			va = compute_va(bird_list, N_birds);
			va_avg += va;
			va_sq_avg += va*va;
		}

		va_avg /= N_measures;
		va_sq_avg /= N_measures;

		fout << eta << '\t' << va_avg << '\t' << std::sqrt(va_sq_avg - va_avg*va_avg) << '\n';

		eta -= .1;
		for (it = bird_list.begin(); it != bird_list.end(); it++) it->set_eta(eta);
	}

	fout.close();

	return 0;
}
