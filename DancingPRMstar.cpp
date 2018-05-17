/*********************************************************************
  * Software License Agreement (BSD License)
  *
  *  Copyright (c) 2013, Willow Garage
  *  All rights reserved.
  *
  *  Redistribution and use in source and binary forms, with or without
  *  modification, are permitted provided that the following conditions
  *  are met:
  *
  *   * Redistributions of source code must retain the above copyright
  *     notice, this list of conditions and the following disclaimer.
  *   * Redistributions in binary form must reproduce the above
  *     copyright notice, this list of conditions and the following
  *     disclaimer in the documentation and/or other materials provided
  *     with the distribution.
  *   * Neither the name of Willow Garage nor the names of its
  *     contributors may be used to endorse or promote products derived
  *     from this software without specific prior written permission.
  *
  *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
  *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
  *  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
  *  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
  *  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
  *  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
  *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
  *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
  *  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
  *  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
  *  POSSIBILITY OF SUCH DAMAGE.
  *********************************************************************/

/* Modified by : A.I Dyk(Donghyuk Kim)*/

#include <omp.h>

#include "ompl/geometric/planners/prm/DancingPRMstar.h"
#include "ompl/base/objectives/PathLengthOptimizationObjective.h"
#include "ompl/base/goals/GoalSampleableRegion.h"
#include "ompl/util/GeometricEquations.h"
#include "ompl/geometric/planners/prm/ConnectionStrategy.h"
#include "ompl/tools/config/SelfConfig.h"
#include "ompl/tools/debug/Profiler.h" // To enable and disable Profiling, please refer to the header file.
#include "ompl/base/spaces/RealVectorStateSpace.h"

#include <boost/lambda/bind.hpp>
#include <boost/graph/astar_search.hpp>
#include <boost/graph/lookup_edge.hpp>
#include <boost/foreach.hpp>
#include <fstream>
#include <queue>
#include <stack>
#include <utility>
#include <numeric>

#define foreach BOOST_FOREACH

#define PI 3.14159265358979

namespace ompl {
namespace magic {
/** \brief The number of nearest neighbors to consider by
      default in the construction of the PRM roadmap */
static const unsigned int DEFAULT_NEAREST_NEIGHBORS_LAZY = 5;

/** \brief When optimizing solutions with lazy planners, this is the minimum
      number of path segments to add before attempting a new optimized solution
          extraction */
static const unsigned int MIN_ADDED_SEGMENTS_FOR_LAZY_OPTIMIZATION = 5;
}
}

ompl::geometric::DancingPRMstar::DistanceField::DistanceField(const Vertex &from, const Vertex &to,
                                                          DancingPRMstar *parent, double resolution)
  : parent_(parent), resolution_(resolution) {
  base::State *state_from = parent->stateProperty_[from], *state_to = parent->stateProperty_[to];

  double cardDbl = static_cast<double>(parent->nn_->size() + 1u);
  double r_nn = std::min(parent_->maxDistance_, parent_->r_rrgConstant_ * std::pow(log(cardDbl) / cardDbl,
                                                                 1 / static_cast<double>(parent_->si_->getStateDimension())));

  std::vector<Vertex> locals;
  locals.push_back(from);
  locals.push_back(to);

  for (auto local = locals.begin(); local != locals.end(); ++local) {
    std::vector<DancingPRMstar::type_neighbor> &neighbors = parent->neighborProperty_[*local];
    for (uint i = 0; i < neighbors.size();) {
      if (neighbors[i].second > parent->max_length_optimized_
        || neighbors[i].second > 2 * r_nn) {
        std::swap(neighbors[i], neighbors.back());
        neighbors.pop_back();
        continue;
      }
      neighbor_list_.insert(neighbors[i].first);
      i++;
    }
  }

  neighbor_list_.insert(from);
  neighbor_list_.insert(to);
}

ompl::geometric::DancingPRMstar::DistanceField::~DistanceField() {
  hash_map_.clear();
}

ompl::geometric::DancingPRMstar::DistanceFieldNode
ompl::geometric::DancingPRMstar::DistanceField::computeDistance(const Eigen::VectorXd st) {
  base::State *state_from = parent_->si_->allocState();

  std::vector<double> v_from(st.data(), st.data() + st.size());
  parent_->si_->getStateSpace()->copyFromReals(state_from, v_from);

  bool isInside = false;
  double closest_dist = std::numeric_limits<double>::infinity();
  Vertex closest_neighbor = nullptr;
  Eigen::VectorXd closest_gradient_vector(Eigen::VectorXd::Zero(v_from.size()));
  Eigen::VectorXd gradient_vector(Eigen::VectorXd::Zero(v_from.size()));

  for (auto it = neighbor_list_.begin(); it != neighbor_list_.end(); ++it) {
    Vertex neighbor = *it;
    double dist_to_neighbor = parent_->distanceFunction(
          parent_->stateProperty_[neighbor], state_from);
    double compensated_distance = dist_to_neighbor - parent_->compensatedRadius(neighbor);
    isInside |= (dist_to_neighbor < 0.0);
    if (isInside) break;

    if (closest_dist > compensated_distance) {
      std::vector<double> v_to;
      parent_->si_->getStateSpace()->copyToReals(v_to, parent_->stateProperty_[neighbor]);

      for (uint i = 0; i < v_from.size(); i++) {
        gradient_vector(i) = (v_to[i] - v_from[i]) / dist_to_neighbor;
      }
      closest_dist = compensated_distance;
      closest_neighbor = neighbor;
    }
  }
  parent_->si_->freeState(state_from);

  if (isInside || closest_neighbor == nullptr) {
    DistanceFieldNode dfn(0.0, Eigen::VectorXd());
    return dfn;
  } else {
    std::vector<double> v_to;
    parent_->si_->getStateSpace()->copyToReals(v_to, parent_->stateProperty_[closest_neighbor]);

    for (uint i = 0; i < v_from.size(); i++) {
      gradient_vector(i) = (v_to[i] - v_from[i]) / (closest_dist + parent_->compensatedRadius((closest_neighbor)));
    }

    DistanceFieldNode dfn(closest_dist + 0.01, gradient_vector);
    return dfn;
  }
}

ompl::geometric::DancingPRMstar::DancingPRMstar(const base::SpaceInformationPtr &si,
    bool starStrategy) :
  base::Planner(si, "DancingPRMstar"),
  starStrategy_(starStrategy),
  userSetConnectionStrategy_(false),
  maxDistance_(0.0),
  indexProperty_(boost::get(boost::vertex_index_t(), g_)),
  stateProperty_(boost::get(vertex_state_t(), g_)),
  radiusProperty_(boost::get(vertex_radius_t(), g_)),
  witnessProperty_(boost::get(vertex_witness_t(), g_)),
  costProperty_(boost::get(vertex_cost_t(), g_)),
  childrenProperty_(boost::get(vertex_children_t(), g_)),
  predecessorProperty_(boost::get(boost::vertex_predecessor_t(), g_)),
  colorProperty_(boost::get(boost::vertex_color_t(), g_)),
  weightProperty_(boost::get(boost::edge_weight_t(), g_)),
  vertexValidityProperty_(boost::get(vertex_flags_t(), g_)),
  edgeValidityProperty_(boost::get(edge_flags_t(), g_)),
  bestCost_(std::numeric_limits<double>::quiet_NaN()),
  iterations_(0),
  increaseIterations_(0),
  dt_(1.0),
  number_of_attempt_opt_(0),
  number_of_success_opt_(0),
  BisectionCC_(true),
  CFreeSpaceApprox_(true),
  CObstacleApprox_(true),
  BiInheritance_(true),
  DancingOpt_(false),
  TestOpt_(false),
  number_of_configuration_(10u),
  rewireFactor_(1.1),
  max_distance_field_value_(0.1),
  gain_(20.0),
  compensationFactor_(0.0),
  projectionFactor_(0.0) {
  specs_.recognizedGoal = base::GOAL_SAMPLEABLE_REGION;
  specs_.approximateSolutions = false;
  specs_.optimizingPaths = true;

  Planner::declareParam<bool>("BisectionCC", this,
                              &DancingPRMstar::setBisectionCC, std::string("."));
  Planner::declareParam<bool>("CFreeSpaceApprox", this,
                              &DancingPRMstar::setCFreeSpaceApprox, std::string("."));
  Planner::declareParam<bool>("CObstacleApprox", this,
                              &DancingPRMstar::setCObstacleApprox, std::string("."));
  Planner::declareParam<bool>("BiInheritance", this,
                              &DancingPRMstar::setBiInheritance, std::string("."));
  Planner::declareParam<bool>("DancingOpt", this,
                              &DancingPRMstar::setDancingOpt, std::string("."));
  Planner::declareParam<bool>("TestOpt", this,
                              &DancingPRMstar::setTestOpt, std::string("."));
  Planner::declareParam<unsigned int>("LengthOfOptTrajectory", this,
                              &DancingPRMstar::setLengthOfOptTrajectory, std::string("."));
  Planner::declareParam<double>("RewireFactor", this,
                              &DancingPRMstar::setRewireFactor, std::string("."));
  Planner::declareParam<double>("MaxDistanceFieldValue", this,
                              &DancingPRMstar::setMaxDistanceFieldValue, std::string("."));
  Planner::declareParam<double>("MaxLengthOptimized", this,
                              &DancingPRMstar::setMaxLengthOptimized, std::string("."));
  Planner::declareParam<double>("MaxMagnitudeOfOptimizedTrajectory", this,
                              &DancingPRMstar::setMaxMagnitudeOfOptimizedTrajectory, std::string("."));
  Planner::declareParam<double>("Gain", this,
                              &DancingPRMstar::setGain, std::string("."));
  Planner::declareParam<double>("CompensationFactor", this,
                              &DancingPRMstar::setCompensationFactor, std::string("."));
  Planner::declareParam<double>("ProjectionFactor", this,
                              &DancingPRMstar::setProjectionFactor, std::string("."));
  Planner::declareParam<double>("epsilon", this,
                              &DancingPRMstar::setEpsilon, std::string("."));
  Planner::declareParam<unsigned int>("MaxOptimizationIteration", this,
                              &DancingPRMstar::setMaxOptimizationIteration, std::string("."));

  if (!starStrategy_) {
    Planner::declareParam<unsigned int>("max_nearest_neighbors", this,
                                       &DancingPRMstar::setMaxNearestNeighbors, std::string("8:1000"));
  }

  addPlannerProgressProperty("iterations INTEGER",
                             boost::bind(&DancingPRMstar::getIterationCount, this));
  addPlannerProgressProperty("best cost REAL",
                             boost::bind(&DancingPRMstar::getBestCost, this));
  addPlannerProgressProperty("milestone count INTEGER",
                             boost::bind(&DancingPRMstar::getMilestoneCountString, this));
  addPlannerProgressProperty("edge count INTEGER",
                             boost::bind(&DancingPRMstar::getEdgeCountString, this));
}

ompl::geometric::DancingPRMstar::~DancingPRMstar() {
}

void ompl::geometric::DancingPRMstar::setup() {
  Planner::setup();
  tools::SelfConfig sc(si_, getName());
  sc.configurePlannerRange(maxDistance_);

  if (!nn_) {
    nn_.reset(tools::SelfConfig::getDefaultNearestNeighbors<Vertex>(this));
    // nn_.reset(new ompl::NearestNeighborsGNAT<Vertex>());
    nn_->setDistanceFunction(boost::bind(&DancingPRMstar::distanceFunctionVertex, this, _1, _2));
  }

  if (!connectionStrategy_) {
    if (starStrategy_) {
      connectionStrategy_ = KStarStrategy<Vertex>(boost::bind(
                              &DancingPRMstar::milestoneCount, this), nn_, si_->getStateDimension());
    } else {
      connectionStrategy_ = KBoundedStrategy<Vertex>
                            (magic::DEFAULT_NEAREST_NEIGHBORS_LAZY, maxDistance_, nn_);
    }
  }

  double dimDbl = static_cast<double>(si_->getStateDimension());
  double prunedMeasure_ = si_->getSpaceMeasure();

  k_rrgConstant_ = rewireFactor_ * boost::math::constants::e<double>() + (boost::math::constants::e<double>() /
                  (double)si_->getStateDimension());

  // r_rrg > 2*(1+1/d)^(1/d)*(measure/ballvolume)^(1/d)
  r_rrgConstant_ = rewireFactor_ * 2.0 * std::pow((1.0 + 1.0/dimDbl) * (prunedMeasure_ / unitNBallMeasure(si_->getStateDimension())), 1.0 / dimDbl);

  OMPL_INFORM("r_rrgConstant_ : %lf", r_rrgConstant_);

  // Setup optimization objective
  //
  // If no optimization objective was specified, then default to
  // optimizing path length as computed by the distance() function
  // in the state space.
  if (pdef_) {
    if (pdef_->hasOptimizationObjective()) {
      opt_ = pdef_->getOptimizationObjective();
    } else {
      opt_.reset(new base::PathLengthOptimizationObjective(si_));

      if (!starStrategy_) {
        opt_->setCostThreshold(opt_->infiniteCost());
      }
    }
  } else {
    OMPL_INFORM("%s: problem definition is not set, deferring setup completion...",
                getName().c_str());
    setup_ = false;
  }

  sampler_ = si_->allocStateSampler();

  // Compensation Ratio
  compensationRatioCache_.resize(1000001);
  for (uint i = 1; i <= 1000000; i++) {
#ifdef COMPENSATION
    compensationRatioCache_[i] = std::max(0.0, 1.0 - (compensationFactor_ * std::pow(log(i) / static_cast<double>(i),
                                                         1 / static_cast<double>(si_->getStateSpace()->getDimension()))));
#else
    compensationRatioCache_[i] = compensationFactor_;
#endif
  }
}

void ompl::geometric::DancingPRMstar::setRange(double distance) {
  maxDistance_ = distance;

  if (!userSetConnectionStrategy_) {
    connectionStrategy_.clear();
  }

  if (isSetup()) {
    setup();
  }
}

void ompl::geometric::DancingPRMstar::setMaxNearestNeighbors(unsigned int k) {
  if (starStrategy_) {
    throw Exception("Cannot set the maximum nearest neighbors for " + getName());
  }

  if (!nn_) {
    nn_.reset(tools::SelfConfig::getDefaultNearestNeighbors<Vertex>(this));
    nn_->setDistanceFunction(boost::bind(&DancingPRMstar::distanceFunctionVertex, this, _1, _2));
  }

  if (!userSetConnectionStrategy_) {
    connectionStrategy_.clear();
  }

  if (isSetup()) {
    setup();
  }
}

void ompl::geometric::DancingPRMstar::setProblemDefinition(
  const base::ProblemDefinitionPtr &pdef) {
  Planner::setProblemDefinition(pdef);
  clearQuery();
}

void ompl::geometric::DancingPRMstar::clearQuery() {
  startM_.clear();
  goalM_.clear();
  pis_.restart();
}

void ompl::geometric::DancingPRMstar::clear() {
  Planner::clear();
  freeMemory();

  if (nn_) {
    nn_->clear();
  }

  clearQuery();

  iterations_ = 0;
  bestCost_ = base::Cost(std::numeric_limits<double>::quiet_NaN());
}

void ompl::geometric::DancingPRMstar::freeMemory() {
  foreach (Vertex v, boost::vertices(g_)) {
    si_->freeState(stateProperty_[v]);
  }

  g_.clear();
}

// Add newly sampled vertex and its adjancency edges connected to neigh neighbors.
ompl::geometric::DancingPRMstar::Vertex ompl::geometric::DancingPRMstar::addMilestone(
  base::State *state, bool isChecked) {
  Vertex m = boost::add_vertex(g_);
  stateProperty_[m] = state;
  radiusProperty_[m] = std::numeric_limits<double>::infinity();
  costProperty_[m] = std::numeric_limits<double>::infinity();
  childrenProperty_[m] = new std::vector<Vertex>();
  witnessProperty_[m] = nullptr;
  predecessorProperty_[m] = nullptr;
  colorProperty_[m] = 0;
  vertexValidityProperty_[m] = (isChecked) ? VALIDITY_TRUE : VALIDITY_UNKNOWN;

  // Which milestones will we attempt to connect to?
  // const std::vector<Vertex> &neighbors = connectionStrategy_(m);

  std::vector<Vertex> neighbors;
  std::vector<double> neighbors_costs;
  const unsigned int max_number_of_neighbors = std::ceil(k_rrgConstant_ * log(static_cast<double>(milestoneCount()) + 1u));

  neighbors.reserve(max_number_of_neighbors);
  neighbors_costs.reserve(max_number_of_neighbors);

  ompl::tools::Profiler::Begin("Nearest Neighbor Search");
  nn_->nearestK(m, max_number_of_neighbors, neighbors);
  ompl::tools::Profiler::End("Nearest Neighbor Search");

  ompl::tools::Profiler::Begin("Configuration-Free Space Approximation");
  // Inherit witness from near neighbors.
  foreach (Vertex n, neighbors) {
    if (witnessProperty_[n] != nullptr) {
      double dist = distanceFunction(witnessProperty_[n], stateProperty_[m]);
      // No directional information for radiusProperty...
      if (dist < radiusProperty_[m]) {
        witnessProperty_[m] = witnessProperty_[n];
        radiusProperty_[m] = dist;
      }
    }
  }
  ompl::tools::Profiler::End("Configuration-Free Space Approximation");

  // Witness Propagation
  foreach (Vertex n, neighbors) {
    ompl::tools::Profiler::Begin("Configuration-Free Space Approximation");
    // Symmetricity assumed here, can't use si_->distance here.
    const double weight = distanceFunction(stateProperty_[n], stateProperty_[m]);
    ompl::base::Cost cost_weight(weight);
    const Graph::edge_property_type properties(cost_weight);

    neighborProperty_[m].push_back(type_neighbor(n, weight));
    neighborProperty_[n].push_back(type_neighbor(m, weight));

    ompl::tools::Profiler::End("Configuration-Free Space Approximation");

    const Edge &e = boost::add_edge(n, m, properties, g_).first;
    edgeValidityProperty_[e] = VALIDITY_UNKNOWN;
    edgeTrajectoryProperty_[e] = nullptr;
  }

  if (witnessProperty_[m] != nullptr) {
    ompl::tools::Profiler::Begin("Configuration-Free Space Approximation");
    foreach (Vertex n, neighbors) {
      double dist = distanceFunction(witnessProperty_[m], stateProperty_[n]);

      if (dist < radiusProperty_[n]) {
        witnessProperty_[n] = witnessProperty_[m];
        radiusProperty_[n] = dist;
      }
    }
    ompl::tools::Profiler::End("Configuration-Free Space Approximation");
  }
  nn_->add(m);

  return m;
}

ompl::base::PlannerStatus ompl::geometric::DancingPRMstar::solve(const
                                                             base::PlannerTerminationCondition &ptc) {
  // Initial checkup for start/goal configurations.
  checkValidity();

  number_of_collision_checks_ = 0;

  // Add the valid start states as milestones
  while (const base::State *st = pis_.nextStart()) {
    Vertex st_vert = addMilestone(si_->cloneState(st));
    costProperty_[st_vert] = 0.0; // Initialize with 0 cost.
    startM_.push_back(st_vert);
  }

  if (startM_.size() == 0) {
    OMPL_ERROR("%s: There are no valid initial states!", getName().c_str());
    return base::PlannerStatus::INVALID_START;
  }

  base::GoalSampleableRegion *goal = dynamic_cast<base::GoalSampleableRegion *>
      (pdef_->getGoal().get());

  if (!goal) {
    OMPL_ERROR("%s: Unknown type of goal", getName().c_str());
    return base::PlannerStatus::UNRECOGNIZED_GOAL_TYPE;
  }

  if (!goal->couldSample()) {
    OMPL_ERROR("%s: Insufficient states in sampleable goal region",
               getName().c_str());
    return base::PlannerStatus::INVALID_GOAL;
  }

  // Ensure there is at least one valid goal state
  if (goal->maxSampleCount() > goalM_.size() || goalM_.empty()) {
    const base::State *st = goalM_.empty() ? pis_.nextGoal(ptc) : pis_.nextGoal();

    if (st) {
      goalM_.push_back(addMilestone(si_->cloneState(st)));
    }

    if (goalM_.empty()) {
      OMPL_ERROR("%s: Unable to find any valid goal states", getName().c_str());
      return base::PlannerStatus::INVALID_GOAL;
    }
  }

  unsigned long int nrStartStates = boost::num_vertices(g_);
  OMPL_INFORM("%s: Starting planning with %lu states already in datastructure",
              getName().c_str(), nrStartStates);

  bestCost_ = opt_->infiniteCost();

  ompl::tools::Profiler::Clear();
  ompl::tools::Profiler::Begin("Total");

  base::State *workState = si_->allocState();
  bool fullyOptimized = false;
  bool someSolutionFound = false;
  base::PathPtr bestSolution;

  Vertex x_obs;
  if (CObstacleApprox_) {
    x_obs = boost::add_vertex(g_);
    radiusProperty_[x_obs] = std::numeric_limits<double>::infinity();
  }

  // Grow roadmap in lazy fashion -- add edges without checking validity
  while (ptc == false) {
    ++iterations_;
    sampler_->sampleUniform(workState);

    // A.I : Pruning can be done here prior to checking its validity.
    /*
    if (!isPromising(workState)) {
      continue;
    }*/

    // A.I : Yet wasting the collision information for sample right now...
    ompl::tools::Profiler::Begin("Vertex Collision Checking");
    number_of_collision_checks_ += 1;
    if (!si_->isValid(workState)) {
      ompl::tools::Profiler::End("Vertex Collision Checking");
      if (CObstacleApprox_) {
        stateProperty_[x_obs] = si_->cloneState(workState);
        Vertex x_closest = nn_->nearest(x_obs);

        ompl::tools::Profiler::Begin("Configuration-Free Space Approximation");
        std::vector<DancingPRMstar::type_neighbor> &neighbors = neighborProperty_[x_closest];

        for (int i = 0; i < (int)neighbors.size(); i++) {
          Vertex neighbor_vertex = neighbors[i].first;
          double dist = distanceFunction(stateProperty_[neighbor_vertex], stateProperty_[x_obs]);

          if (dist < radiusProperty_[neighbor_vertex]) {
            radiusProperty_[neighbor_vertex] = dist;
            witnessProperty_[neighbor_vertex] = stateProperty_[x_obs];
          }
        }
        ompl::tools::Profiler::End("Configuration-Free Space Approximation");
      }
      continue;
    } else {
      ompl::tools::Profiler::End("Vertex Collision Checking");
    }

    // Add collision-free vertices.
    Vertex addedVertex = addMilestone(si_->cloneState(workState), true);

    // DSPT update.
    Decrease(addedVertex);

    // Only support a single pair of start and goal node.
    Vertex startV = startM_[0];
    Vertex goalV = goalM_[0];
    base::PathPtr solution;

    do {
      if (predecessorProperty_[goalV] == nullptr || bestCost_.value() <= costProperty_[goalV])
        break;
      solution = constructSolution(startV, goalV);
    } while (!solution);

    if (solution) {
      someSolutionFound = true;
      base::Cost c(costProperty_[goalV]);

      // A.I : If it is a non-optimal planner, it will be terminated at here,
      // will keep iterating otherwise.
      if (opt_->isSatisfied(c)) {
        fullyOptimized = true;
        bestSolution = solution;
        bestCost_ = c;
        break;
      } else {
        if (opt_->isCostBetterThan(c, bestCost_)) {
          bestSolution = solution;
          OMPL_INFORM("%.6lf -> %.6lf", bestCost_.value(), c.value());
          bestCost_ = c;
        }
      }
    }
  }
  si_->freeState(workState);

  if (bestSolution) {
    base::PlannerSolution psol(bestSolution);
    psol.setPlannerName(getName());
    // If the solution was optimized, we mark it as such
    psol.setOptimized(opt_, bestCost_, fullyOptimized);
    pdef_->addSolutionPath(psol);
  }

  ompl::tools::Profiler::End("Total");

  if (CObstacleApprox_) boost::remove_vertex(x_obs, g_);

  OMPL_INFORM("%lf", bestCost_.value());

  double avg_degree = 0.0;
  BGL_FORALL_VERTICES(vert, g_, Graph)
      avg_degree += boost::out_degree(vert, g_);
  avg_degree /= boost::num_vertices(g_);

  OMPL_INFORM("%s: Created %u vertices and %u edges with avg. degree %lf.", getName().c_str(),
              boost::num_vertices(g_) - 1, boost::num_edges(g_), avg_degree);

  OMPL_INFORM("Opt success rate : %.2lf \% (%u/%u)", (double)number_of_success_opt_ /
             number_of_attempt_opt_ * 100.0, number_of_success_opt_, number_of_attempt_opt_);

  return bestSolution ? base::PlannerStatus::EXACT_SOLUTION :
                        base::PlannerStatus::TIMEOUT;
}

// Check optimized trajectory, it can over check regardless of collision checking resolution.
double ompl::geometric::DancingPRMstar::checkTrajectory(const std::vector<base::State*> &v, const Vertex &v1, const Vertex &v2,
                                                    const double weight) {
  ompl::tools::Profiler::ScopedBlock _profiler("Edge Collision Checking");
  double trajectory_cost = 0.0f;
  auto prev = v.begin();
  bool witness_update = false;
  std::vector<base::State*> temp;
  trajectory_cost += distanceFunction(stateProperty_[v1], v.front());
  temp.push_back(stateProperty_[v1]);
  temp.push_back(*prev);
  for (auto it = prev + 1; it != v.end(); ++it) {
    trajectory_cost += distanceFunction(*prev, *it);
    temp.push_back(*it);
    prev = it;
  }
  temp.push_back(stateProperty_[v2]);
  trajectory_cost += distanceFunction(v.back(), stateProperty_[v2]);

  if (trajectory_cost > max_magnitude_of_optimized_trajectory_ * weight) return -1.0;

  // Compute validSegmentCount for arbitrary length of trajectory.
  /*
  unsigned int nd = si_->getStateSpace()->getValidSegmentCountFactor() *
                trajectory_cost / si_->getStateSpace()->getLongestValidSegmentLength();
  double measurable_length = trajectory_cost / static_cast<double>(nd);
  */
  double measurable_length = si_->getStateSpace()->getLongestValidSegmentLength() /
      si_->getStateSpace()->getValidSegmentCountFactor();
  double accumulated_cost = measurable_length;

  base::State *pivot = si_->allocState();
  prev = temp.begin();
  auto it = prev + 1;

  while (it != temp.end()) {
    double dist = distanceFunction(*prev, *it);

    while (accumulated_cost < dist) {
      si_->getStateSpace()->interpolate(*prev, *it, accumulated_cost / dist, pivot);

      number_of_collision_checks_++;
      if (si_->isValid(pivot)) {
        // Collision test passed
      } else {
        if (CFreeSpaceApprox_) {
          ompl::tools::Profiler::Begin("Configuration-Free Space Approximation");
          dist = distanceFunction(pivot, stateProperty_[v1]);

          // A.I : Update witness and radius.
          if (radiusProperty_[v1] > dist) {
            witnessProperty_[v1] = pivot;
            radiusProperty_[v1] = dist;
          }
          dist = distanceFunction(pivot, stateProperty_[v2]);

          if (radiusProperty_[v2] > dist) {
            witnessProperty_[v2] = pivot;
            radiusProperty_[v2] = dist;
          }

          std::vector<DancingPRMstar::type_neighbor> &neighbors_v1 = neighborProperty_[v1];
          for (int i = 0; i < (int)neighbors_v1.size(); i++) {
            Vertex neighbor_vertex = neighbors_v1[i].first;
            double dist2 = distanceFunction(pivot, stateProperty_[neighbor_vertex]);

            if (radiusProperty_[neighbor_vertex] > dist2) {
              witnessProperty_[neighbor_vertex] = pivot;
              radiusProperty_[neighbor_vertex] = dist2;
            }
          }

          std::vector<DancingPRMstar::type_neighbor> &neighbors_v2 = neighborProperty_[v2];
          for (int i = 0; i < (int)neighbors_v2.size(); i++) {
            Vertex neighbor_vertex = neighbors_v2[i].first;
            double dist2 = distanceFunction(pivot, stateProperty_[neighbor_vertex]);

            if (radiusProperty_[neighbor_vertex] > dist2) {
              witnessProperty_[neighbor_vertex] = pivot;
              radiusProperty_[neighbor_vertex] = dist2;
            }
          }

          ompl::tools::Profiler::End("Configuration-Free Space Approximation");
        }
        return -1.0;
      }

      accumulated_cost += measurable_length;
    }
    accumulated_cost -= dist;

    prev = it;
    ++it;
  }
  // if (!witness_update)
  si_->freeState(pivot);

  return trajectory_cost;
}

// A.I : Check adaptive edge collision checking from v1 to v2 and return witness closest from v1.
bool ompl::geometric::DancingPRMstar::checkMotionAdaptively(const Vertex &v1, const Vertex &v2,
                                                         const double weight, std::vector<base::State*> *optimized) {
  /* Assume motion starts/ends in a valid configuration so v1/v2 are valid */

  const double r1 = radiusProperty_[v1], r2 = radiusProperty_[v2];
  bool result = true;

  if (result) {
    // valid_++;
  } else {
    // invalid_++;
  }

  return result;
}

// A.I : Check edge collision checking from v1 to v2 and update witness.
bool ompl::geometric::DancingPRMstar::checkMotion(const Vertex &v1, const Vertex &v2,
    const double weight) {
  ompl::tools::Profiler::ScopedBlock _profiler("Edge Collision Checking");
  /* Assume motion starts/ends in a valid configuration so v1/v2 are valid */

  bool result = true, witness_update = false;
  const base::State *s1 = stateProperty_[v1], *s2 = stateProperty_[v2];

  int nd = si_->getStateSpace()->validSegmentCount(s1, s2);
  double from_radius = (double)nd * radiusProperty_[v1] / weight;
  double to_radius = (double)nd * (1.0 - radiusProperty_[v2] / weight);

  if (nd > 1) {
    base::State *pivot = si_->allocState();
    /* Temporary storage for the checked state */
    if (!BisectionCC_) {
      // At j == 0 and nd -1 is corresponds to s1 and s2, respectively.
      for (int j = 1 ; j < nd ; j++) {
        si_->getStateSpace()->interpolate(s1, s2, (double)j / (double)nd, pivot);

        number_of_collision_checks_++;
        if (!si_->isValid(pivot)) {
          if (CFreeSpaceApprox_) {
            ompl::tools::Profiler::Begin("Configuration-Free Space Approximation");
            // A.I : Update witness and radius.
            // if ((double)j <= from_radius)  {
            if (radiusProperty_[v1] > ((double)j / (double)nd) * weight) {
              witnessProperty_[v1] = pivot;
              radiusProperty_[v1] = ((double)j / (double)nd) * weight;
            }

            // if (to_radius <= (double)j) {
            if (radiusProperty_[v2] > (1.0 - (double)j / (double)nd) * weight) {
              witnessProperty_[v2] = pivot;
              radiusProperty_[v2] = (1.0 - (double)j / (double)nd) * weight;
            }

            std::vector<DancingPRMstar::type_neighbor> &neighbors_v1 = neighborProperty_[v1];
            for (int i = 0; i < (int)neighbors_v1.size(); i++) {
              Vertex neighbor_vertex = neighbors_v1[i].first;
              double dist2 = distanceFunction(pivot, stateProperty_[neighbor_vertex]);

              if (radiusProperty_[neighbor_vertex] > dist2) {
                witnessProperty_[neighbor_vertex] = pivot;
                radiusProperty_[neighbor_vertex] = dist2;
              }
            }

            std::vector<DancingPRMstar::type_neighbor> &neighbors_v2 = neighborProperty_[v2];
            for (int i = 0; i < (int)neighbors_v2.size(); i++) {
              Vertex neighbor_vertex = neighbors_v2[i].first;
              double dist2 = distanceFunction(pivot, stateProperty_[neighbor_vertex]);

              if (radiusProperty_[neighbor_vertex] > dist2) {
                witnessProperty_[neighbor_vertex] = pivot;
                radiusProperty_[neighbor_vertex] = dist2;
              }
            }
            ompl::tools::Profiler::End("Configuration-Free Space Approximation");
          }

          return false;
        }
      }
    } else {
      std::queue<std::pair<unsigned int, unsigned int> > q;
      q.push(std::make_pair(1, nd - 1));

      while (!q.empty()) {
        auto range = q.front();
        unsigned int mid;

        mid = (range.first + range.second) / 2;
        si_->getStateSpace()->interpolate(s1, s2, (double)mid / (double)nd, pivot);

        number_of_collision_checks_++;
        if (!si_->isValid(pivot)) {
          if (CFreeSpaceApprox_) {
            ompl::tools::Profiler::Begin("Configuration-Free Space Approximation");
            int is_outside = 0;
            // A.I : Update witness and radius.
            double distance_to_obs = ((double)mid / (double)nd) * weight; // from v1
            if (radiusProperty_[v1] > distance_to_obs) {
              witnessProperty_[v1] = pivot;
              radiusProperty_[v1] = distance_to_obs;
              is_outside += 1;
            }

            distance_to_obs = (1.0 - (double)mid / (double)nd) * weight; // from v2
            if (radiusProperty_[v2] > distance_to_obs) {
              witnessProperty_[v2] = pivot;
              radiusProperty_[v2] = distance_to_obs;
              is_outside += 1;
            }

            std::vector<DancingPRMstar::type_neighbor> &neighbors_v1 = neighborProperty_[v1];
            for (int i = 0; i < (int)neighbors_v1.size(); i++) {
              Vertex neighbor_vertex = neighbors_v1[i].first;
              double dist2 = distanceFunction(pivot, stateProperty_[neighbor_vertex]);

              if (radiusProperty_[neighbor_vertex] > dist2) {
                witnessProperty_[neighbor_vertex] = pivot;
                radiusProperty_[neighbor_vertex] = dist2;
              }
            }

            std::vector<DancingPRMstar::type_neighbor> &neighbors_v2 = neighborProperty_[v2];
            for (int i = 0; i < (int)neighbors_v2.size(); i++) {
              Vertex neighbor_vertex = neighbors_v2[i].first;
              double dist2 = distanceFunction(pivot, stateProperty_[neighbor_vertex]);

              if (radiusProperty_[neighbor_vertex] > dist2) {
                witnessProperty_[neighbor_vertex] = pivot;
                radiusProperty_[neighbor_vertex] = dist2;
              }
            }
            ompl::tools::Profiler::End("Configuration-Free Space Approximation");
          }

          return false;
        }

        q.pop();
        if (range.first < mid)
          q.push(std::make_pair(range.first, mid - 1));
        if (mid < range.second) {
          q.push(std::make_pair(mid + 1, range.second));
        } // if mid == first, no more recursion.
      }
    }

    si_->freeState(pivot);  // If it is collision-free, no witness.
  }

  return result;
}

bool ompl::geometric::DancingPRMstar::isPromising(const base::State *s) {
  if (!opt_->isFinite(bestCost_))
    return true;

  double dist = si_->distance(stateProperty_[startM_[0]], s) +
                si_->distance(s, stateProperty_[goalM_[0]]);
  return dist < bestCost_.value();
}

// outedge, inedge? - doesn't matter, need to scan all the neighbors.
void ompl::geometric::DancingPRMstar::Decrease(const Vertex &v) {
  ompl::tools::Profiler::ScopedBlock _profiler("Decrease");
  typedef std::pair<double, Vertex> weight_vertex;
  std::priority_queue<weight_vertex> pq;

  // Initialize cost of v, i.e., finding best parent vertex in G(g_).
  BGL_FORALL_OUTEDGES(v, e, g_, Graph) {
    Vertex w = target(e, g_);
    double weight = weightProperty_[e].value();

    if (costProperty_[v] > costProperty_[w] + weight) {
      predecessorProperty_[v] = w;
      costProperty_[v] = costProperty_[w] + weight;
    }
  }

  // No need to invoke cancelAdoption since v is newly sampled.
  if (predecessorProperty_[v] != nullptr) {
    childrenProperty_[predecessorProperty_[v]]->push_back(v);
  }

  // At this point, v has a best parent. From now on construct its subtree of descendants.

  pq.push(weight_vertex(-costProperty_[v], v)); // Invert the cost value for mimicking min-heap.

  while (!pq.empty()) {
    weight_vertex top = pq.top();
    pq.pop();

    double cost = -top.first; // Invert the cost value to be like min-heap.
    Vertex vert = top.second;

    if (cost > costProperty_[vert]) {
      continue;
    }

    BGL_FORALL_OUTEDGES(vert, e, g_, Graph) {
      Vertex w = target(e, g_);
      double weight = weightProperty_[e].value();
      double cost_w = costProperty_[w];

      if (cost_w > cost + weight) {
        costProperty_[w] = cost + weight;
        cancelAdoption(w);

        predecessorProperty_[w] = vert;
        childrenProperty_[vert]->push_back(w);
        pq.push(weight_vertex(-costProperty_[w], w));
      }
    }
  }

  // Now, DSPT is stable.
}

#define RED (increaseIterations_ + 1) // I know, it's bad #define.

void ompl::geometric::DancingPRMstar::Increase(const Vertex vs) {
  ompl::tools::Profiler::ScopedBlock _profiler("Increase");
  // <Step 1. Preparation.
  // white is used for color of each vertex without initialization.
  // For each iteration, white is increased by 1, thus we can use it as
  // equal to or less than 'white' means 'white' color, 'red' otherwise.
  increaseIterations_ += 1;

  std::vector<Vertex> reds;
  typedef std::pair<double, Vertex> weight_vertex;
  std::priority_queue<weight_vertex> pq; //  Max-heap by default.

  pq.push(weight_vertex(-costProperty_[vs], vs));  // It works as if it is min-heap.

  // <Step 2. Coloring
  while (!pq.empty()) {
    weight_vertex top = pq.top();
    pq.pop();

    double cost = -top.first;
    Vertex vert = top.second;

    if (cost > costProperty_[vert]) {
      continue;  // Instead of heap-improve
    }

    // Probability to get a pink node? almost impossible.
    /*
    // If there exist a non-red neighbor q of z such that Dist(q) + w_(q, z) = D(z)
    // set pink, that means it can keep the current cost, thus it is not necessary to
    // iterate its children.
    // Otherwise, set red and enqueue all the children of z.
    bool pink_flag = false;
    BGL_FORALL_OUTEDGES(vert, e, g_, Graph) {
      Vertex w = target(e, g_);
      double weight = weightProperty_[e].value();

      if (colorProperty_[w] != RED && costProperty_[w] + weight == cost) {
        // Actually, '<' should not be happened all the time.
        // And even '==' would very rarely occur, but possible.
        pink_flag = true;

        cancelAdoption(vert);
        predecessorProperty_[vert] = w;
        childrenProperty_[w]->push_back(vert);
        break; // If there exsits, take anyone among them.
      }
    }

    if (pink_flag) {
      continue;
    }*/

    colorProperty_[vert] = RED; // Set to 'red'
    reds.push_back(vert);
    // Even with multiple starting red nodes, there will be no re-visit since each starting node is
    // a root node of sub'tree' in DSPT. That is, if statement within for loop might be useless.
    // Someone would be curious, e.g., then why do we have to use priority queue in step2 ?
    // Just for 'pink' case. Yeap. We need to identify all the other parent candidates are red or not
    // prior to checking current node.
    std::vector<Vertex> *children = childrenProperty_[vert];

    for (unsigned int i = 0; i < children->size(); i++) {
      pq.push(weight_vertex(-costProperty_[(*children)[i]], (*children)[i]));
    }
  }

  // 'pq' is empty at here.

  // <Step 3-a. Find best non-red parent for each red node.
  for (unsigned int i = 0; i < reds.size(); i++) {
    // TODO : need to be verified
    // Cost/predecessor initialization.
    costProperty_[reds[i]] = std::numeric_limits<double>::infinity();
    cancelAdoption(reds[i]);

    BGL_FORALL_OUTEDGES(reds[i], e, g_, Graph) {
      Vertex w = target(e, g_);
      double weight = weightProperty_[e].value();

      if (colorProperty_[w] == RED) {
        continue;  // If red, put aside for a while.
      }

      if (costProperty_[reds[i]] > costProperty_[w] + weight) {
        costProperty_[reds[i]] = costProperty_[w] + weight;
        predecessorProperty_[reds[i]] = w;
      }
    }

    if (predecessorProperty_[reds[i]] != nullptr) {
      childrenProperty_[predecessorProperty_[reds[i]]]->push_back(reds[i]);
    }

    pq.push(weight_vertex(-costProperty_[reds[i]], reds[i]));
  }

  // <Step 3-b. Propagate the changes; rewiring for 'red' nodes whether it can replace
  // existing parent node of near neighbors.
  while (!pq.empty()) {
    weight_vertex top = pq.top();
    pq.pop();

    double cost = -top.first;
    Vertex vert = top.second;

    if (costProperty_[vert] < cost) {
      continue;
    }

    BGL_FORALL_OUTEDGES(vert, e, g_, Graph) {
      Vertex w = target(e, g_);
      double weight = weightProperty_[e].value();

      if (colorProperty_[w] != RED) {
        continue;  // If not red, then skip.
      }

      if (cost + weight < costProperty_[w]) {
        costProperty_[w] = cost + weight;

        cancelAdoption(w);

        predecessorProperty_[w] = vert;
        childrenProperty_[vert]->push_back(w);
        pq.push(weight_vertex(-costProperty_[w], w));
      }
    }
  }

  // The end!
  // colorProperty_ is not necessary to be cleansed out, just increase variable, 'RED'.
}

// TODO : sync between children & edges.
void ompl::geometric::DancingPRMstar::cancelAdoption(const Vertex &child) {
  // ompl::tools::Profiler::ScopedBlock _profiler("Remove From Parent");
  if (predecessorProperty_[child] == nullptr)
    return;

  std::vector<Vertex> *children = childrenProperty_[predecessorProperty_[child]];

  for (unsigned int i = 0; i < children->size(); i++) if ((*children)[i] == child) {
    std::swap((*children)[i], children->back());
    children->pop_back();
    break;
  }

  predecessorProperty_[child] = nullptr;
}

// Vertex first.
ompl::base::PathPtr ompl::geometric::DancingPRMstar::constructSolution(
  const Vertex &start, const Vertex &goal) {
  ompl::tools::Profiler::ScopedBlock _profiler("Shortest Path Computation");
  std::vector<Vertex> solution_path;

  // Construct a solution from DSPT.
  for (Vertex vert = goal; vert != nullptr; vert = predecessorProperty_[vert])
    solution_path.push_back(vert);

  if (solution_path.empty())
    return base::PathPtr();

  // Goal = 0, Start = n - 1.
  // TODO : From goal or start ? which one is better?

  auto from = solution_path.rbegin();
  for (auto to = from + 1; to != solution_path.rend(); ++to) {
    Edge e = boost::lookup_edge(*from, *to, g_).first; // Exhaustive search O(E) at worst case.
    unsigned int &evd = edgeValidityProperty_[e];

    if ((evd & VALIDITY_TRUE) == 0) { // Unknown
      bool result = true;
      double weight = weightProperty_[e].value();

      result &= checkMotion(*from, *to, weightProperty_[e].value());

      if (result) {
        evd |= VALIDITY_TRUE;
      } else {
        if (DancingOpt_) {
          double cardDbl = static_cast<double>(nn_->size() + 1u);
          double r_nn = std::min(maxDistance_, r_rrgConstant_ * std::pow(log(cardDbl) / cardDbl,
                                                                         1 / static_cast<double>(si_->getStateDimension())));
          // Optimize when they are close enough
          // to cover the edge by free spheres of neighbors.

          if (2 * r_nn > weight && weight > si_->getStateSpace()->getLongestValidSegmentLength()) {
            std::vector<base::State*> *optimized = optimizeMotion(*from, *to);
            number_of_attempt_opt_ += 1;
            if (optimized != nullptr) {
              double traj_result = checkTrajectory(*optimized, *from, *to, weight);

              if (traj_result > 0) {
                number_of_success_opt_ += 1;
                evd |= VALIDITY_TRUE;
                edgeTrajectoryProperty_[e] = optimized;
                weightProperty_[e] = base::Cost(traj_result);

                for (uint i = 0 ; i < optimized->size(); i++) {
                  trajectory_history_.push_back((*optimized)[i]);
                }
              } else {
                boost::remove_edge(e, g_); // O(log(E/V)) time...
              }

              cancelAdoption(*to);
              Increase(*to);
              return base::PathPtr();
            }
          }
        }
        boost::remove_edge(e, g_); // O(log(E/V)) time...

        cancelAdoption(*to);
        Increase(*to);
        return base::PathPtr();
      }
    }

    from = to;
  }

  PathGeometric *p = new PathGeometric(si_);

  // Feasible path is found, fetch optimized edges if possible.
  for (std::vector<Vertex>::const_reverse_iterator sol = solution_path.rbegin();
       sol != solution_path.rend(); ++sol) {
    p->append(stateProperty_[*sol]);

    if (sol + 1 != solution_path.rend()) {
      Edge e = boost::lookup_edge(*(sol + 1), *(sol), g_).first;
      std::vector<base::State*> *traj = edgeTrajectoryProperty_[e];
      if (traj != nullptr) {
        for (auto it = traj->begin(); it != traj->end(); ++it) {
          p->append(*it);
        }
      }
    }
  }

  return base::PathPtr(p);
}

void ompl::geometric::DancingPRMstar::getPlannerData(base::PlannerData &data)
const {
  Planner::getPlannerData(data);
  if (DancingOpt_) {
    uint i = 0;
    while (i < trajectory_history_.size()) {
      int disguised_float;
      float packet;
      // Iteation through number_of_configuration - 1 for intermediate nodes,
      // + 2 for start/end nodes.
      for (uint j = i; j < i + number_of_configuration_ - 1; j++) {
        data.addEdge(base::PlannerDataVertex(trajectory_history_[j]),
                     base::PlannerDataVertex(trajectory_history_[j + 1]));
        packet = (j == i) ? -1.0 : -2.0;
        memcpy(&disguised_float, &packet, sizeof(int));
        data.tagState(trajectory_history_[j], disguised_float); // Separator.
      }
      packet = -2.0;
      memcpy(&disguised_float, &packet, sizeof(int));
      data.tagState(trajectory_history_[i + number_of_configuration_ - 1], disguised_float);
      i += number_of_configuration_;
    }

    /*
    foreach (const Vertex v, boost::vertices(g_)) {
      // const Vertex v1 = boost::source(e, g_);
      // const Vertex v2 = boost::target(e, g_);
      if (witnessProperty_[v] != nullptr) {
        data.addEdge(base::PlannerDataVertex(stateProperty_[v]),
                     base::PlannerDataVertex(witnessProperty_[v]));

        float radius = static_cast<float>(compensatedRadius(v));
        int disguised_float;
        memcpy(&disguised_float, &radius, sizeof(int));
        data.tagState(stateProperty_[v], disguised_float);

        radius = -3.0;
        memcpy(&disguised_float, &radius, sizeof(int));
        data.tagState(witnessProperty_[v], disguised_float);
      }
    }*/

    /*
    foreach (const Edge e, boost::edges(g_)) {
      const Vertex v1 = boost::source(e, g_);
      const Vertex v2 = boost::target(e, g_);

      if (edgeTrajectoryProperty_[e] == nullptr) { // straight segment
        data.addEdge(base::PlannerDataVertex(stateProperty_[v1]),
                     base::PlannerDataVertex(stateProperty_[v2]));

        float radius = static_cast<float>(compensatedRadius(v1));
        int disguised_float;
        memcpy(&disguised_float, &radius, sizeof(int));
        data.tagState(stateProperty_[v1], disguised_float);

        radius = static_cast<float>(compensatedRadius(v2));
        memcpy(&disguised_float, &radius, sizeof(int));
        data.tagState(witnessProperty_[v2], disguised_float);
      }
    }
    */

    typedef std::pair<Vertex, bool> q_node;
    std::stack<q_node> q;
    q.push(q_node(startM_[0], false));
    while (!q.empty()) {
      q_node top_node = q.top();
      Vertex top = top_node.first;
      q.pop();

      // OMPL_DEBUG("%u:%u", childrenProperty_[top]->size(), predecessorProperty_[top] == nullptr);

      if ((childrenProperty_[top] == nullptr || childrenProperty_[top]->empty() || top_node.second) && predecessorProperty_[top] != nullptr) {
        Edge e = boost::lookup_edge(predecessorProperty_[top], top, g_).first;
        Edge e2 = boost::lookup_edge(top, predecessorProperty_[top], g_).first;
        if (edgeTrajectoryProperty_[e] == nullptr) {
          data.addEdge(base::PlannerDataVertex(stateProperty_[top]),
                    base::PlannerDataVertex(stateProperty_[predecessorProperty_[top]]));

          float radius = static_cast<float>(compensatedRadius(top));
          int disguised_float;
          if (edgeValidityProperty_[e] == VALIDITY_TRUE)
            radius += 10;
          else if (edgeValidityProperty_[e2] == VALIDITY_TRUE)
            radius += 10;
          memcpy(&disguised_float, &radius, sizeof(int));
          data.tagState(stateProperty_[top], disguised_float);

          radius = static_cast<float>(compensatedRadius(predecessorProperty_[top]));
          radius += edgeValidityProperty_[e] == VALIDITY_TRUE ? 10 : 0;
          memcpy(&disguised_float, &radius, sizeof(int));
          data.tagState(witnessProperty_[predecessorProperty_[top]], disguised_float);
        } else {
          std::vector<base::State*> *traj = edgeTrajectoryProperty_[e];
          auto prev = traj->begin();
          for (auto it = prev + 1; it != traj->end(); ++it) {
            data.addEdge(base::PlannerDataVertex(*prev),
                         base::PlannerDataVertex(*it));

            float radius = 100.0f;
            int disguised_float;
            memcpy(&disguised_float, &radius, sizeof(int));
            data.tagState(*prev, disguised_float);
            data.tagState(*it, disguised_float);
            prev = it;
          }
        }
        continue;
      } else if (top_node.second && predecessorProperty_[top] == nullptr) {
        break;
      }
      q.push(q_node(top, true));

      for (uint i = 0; i < childrenProperty_[top]->size(); i++) {
        Vertex next = (*childrenProperty_[top])[i];
        q.push(q_node(next, false));
      }
    }
  }
}

// Local parameters should be merged into DancingPRMstar.
std::vector<ompl::base::State*>* ompl::geometric::DancingPRMstar::optimizeMotion(const Vertex &v1, const Vertex &v2) {
  // The implementation of local trajectory optimizer is not provided due to the
  // patent/project-related issues.

  // For the implementation refer to the
  // CHOMP(Covariant Hamiltonian Optimization for Motion Planning, N. Ratliff et al., ICRA 2009)
  // http://wiki.ros.org/chomp_motion_planner
  // http://wiki.ros.org/stomp_motion_planner
  return nullptr;

  // Belows are HowTo for the use of DistanceField to estimate the distance to the closest
  // approx. C-free space from the given configuration.

  // Initialize local distance_field from the approx. C-free space.
  DistanceField distance_field(v1, v2, this);

  // A query configuration which would be a configuration along the trajectory to be optimized.
  Eigen::VectorXd const xx;

  // During the optimization iteration, we compute f_obs(.) using computeDistance(.).
  DistanceFieldNode dfn = distance_field.computeDistance(xx);

  // Distance to the closest approx. C-free space.
  const double dist_max = max_distance_field_value_;
  const double dist = std::min(dist_max, dfn.dist_);

  // gradient_ for \nabla f_obs computation. It contains a d-dimensional vector toward the closest
  // approx. C-free space.
  Eigen::VectorXd delta = dfn.gradient_;
}

double ompl::geometric::DancingPRMstar::compensatedRadius(const Vertex n) const {
  if (compensationFactor_ == 0.0) return radiusProperty_[n];
  return radiusProperty_[n] * compensationRatioCache_[std::min(1000000ul, (unsigned long)number_of_collision_checks_) + 1u];
}

double ompl::geometric::DancingPRMstar::distanceFunction(const base::State *a, const base::State *b) const {
  return si_->distance(a, b);
}

double ompl::geometric::DancingPRMstar::distanceFunctionVertex(const Vertex a, const Vertex b) const {
  return si_->distance(stateProperty_[a], stateProperty_[b]);
}
