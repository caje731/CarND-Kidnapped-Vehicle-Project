/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h>
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of
	//   x, y, theta and their uncertainties from GPS) and all weights to 1.
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

  default_random_engine generator;
  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);

  num_particles = 70;
  for (int i = 0; i < num_particles; i++) {
    Particle p = Particle();
    p.id = i;
    p.x = dist_x(generator);
    p.y = dist_y(generator);
    p.theta = dist_theta(generator);
    p.weight = 1;
    particles.push_back(p);
    weights.push_back(1);
  }

  is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
  double std_x = std_pos[0];
  double std_y = std_pos[1];
  double std_psi = std_pos[2];

  default_random_engine gen;
  // create a normal (Gaussian) distribution for x, y and psi
  normal_distribution<double> dist_x(0.0, std_x);
  normal_distribution<double> dist_y(0.0, std_y);
  normal_distribution<double> dist_psi(0.0, std_psi);

  if (yaw_rate == 0.0) {
    std::cout << "yaw rate zero" << std::endl;

  } else {
//    std::cout << "velocity " << velocity << std::endl;
//    std::cout << "yaw rate " << yaw_rate << std::endl;

    for (int i = 0; i < num_particles; i++) {
      Particle &p = particles[i];
      double x_d = velocity / yaw_rate * (sin(p.theta + yaw_rate * delta_t) - sin(p.theta));
      double y_d = velocity / yaw_rate * (cos(p.theta) - cos(p.theta + yaw_rate * delta_t));
      double t_d = yaw_rate * delta_t;

      double x_g = dist_x(gen);
      double y_g = dist_y(gen);
      double t_g = dist_psi(gen);

//      std::cout << "particle " << p.id << std::endl;
//      std::cout << "x " << p.x << " + " << x_d << " + " << x_g << std::endl;
//      std::cout << "y " << p.y << " + " << y_d << " + " << y_g << std::endl;
//      std::cout << "t " << p.theta << " + " << t_d << " + " << t_g << std::endl;

      p.x = p.x + x_d + x_g;
      p.y = p.y + y_d + y_g;
      p.theta = p.theta + t_d + t_g;
    }
  }
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
	//   implement this method and use it as a helper during the updateWeights phase.
  for (int i = 0; i < observations.size(); ++i) {
    LandmarkObs &observation = observations[i];

    LandmarkObs closest_landmark;
    double closest_distance = 1000; //Just needs to be bigger than sensor_range!
    for (LandmarkObs map_landmark : predicted) {
      double distance = sqrt(pow(observation.x - map_landmark.x, 2) + pow(observation.y - map_landmark.y, 2));
      if (distance < closest_distance) {
        closest_distance = distance;
        closest_landmark = map_landmark;
      }
    }
    if (closest_distance < 1000) {
      observation.id = closest_landmark.id;
    }
  }
//  cout << "end dataAssociation" << endl;
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html

  // need a position to initialize the coordinate conversion
  for (int i = 0; i < num_particles; i++) {
    Particle &particle = particles[i];
    // measurements in MAP coordinate system (based on particle)
    std::vector<LandmarkObs> observations_map;
    for (LandmarkObs obs : observations) {
      double x_map = particle.x + cos(particle.theta) * obs.x - sin(particle.theta) * obs.y;
      double y_map = particle.y + sin(particle.theta) * obs.x + cos(particle.theta) * obs.y;

      LandmarkObs obs_map = LandmarkObs();
      //obs_map.id = obs.id; id not set yet
      obs_map.x = x_map;
      obs_map.y = y_map;
      observations_map.push_back(obs_map);
    }

    // data association, assume predicted = map and observation = measurement
    // we know we can't sense any landmark beyond sensor_range
    // convert Map single_landmark_s to LandmarkObs?
    std::vector<LandmarkObs> predicted;
    for (Map::single_landmark_s landmark : map_landmarks.landmark_list) {
      double distance = sqrt(pow(particle.x - landmark.x_f, 2) + pow(particle.y - landmark.y_f, 2));
      if (distance > sensor_range) {
//        cout << "Reject " << landmark.id_i << " distance " << distance << endl;
      } else {
        LandmarkObs pred_landmark = LandmarkObs();
        pred_landmark.id = landmark.id_i;
        pred_landmark.x = landmark.x_f;
        pred_landmark.y = landmark.y_f;
        predicted.push_back(pred_landmark);
      }
    }
    // updates id in observations by reference
    dataAssociation(predicted, observations_map);

    // do for each particle
    double prob = 1.0;
    for (LandmarkObs observation : observations_map) {
      if (observation.id == 0) {
        prob *= 0;
      } else {
//        cout << "observation: " << observation.id << endl;
        Map::single_landmark_s landmark = map_landmarks.landmark_list[observation.id - 1];
//        cout << "map: " << landmark.id_i << endl;

        double x_diff = observation.x - landmark.x_f;
        double y_diff = observation.y - landmark.y_f;
        double gaussian = (1.0 / (2 * M_PI * std_landmark[0] * std_landmark[1])) *
                          exp(-0.5 * ((1.0 / pow(std_landmark[0], 2)) * x_diff * x_diff +
                                      (1.0 / pow(std_landmark[1], 2)) * y_diff * y_diff));

//        cout << "gaussian: " << gaussian << endl;
        prob *= gaussian;
      }
    }
//    cout << i << " final prob: " << prob << endl;
    particle.weight = prob;
    weights[i] = prob;
  }
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight.
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

  std::vector<Particle> new_particles;
  std::random_device rd;
  std::mt19937 gen(rd());
  std::discrete_distribution<> d(weights.begin(), weights.end());
  for (int i = 0; i < num_particles; i++) {
    new_particles.push_back(particles[d(gen)]);
  }
  particles = new_particles;

}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
