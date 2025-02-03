# GA Worker Assignment 🚀

Hello and welcome! 👋

**GA Worker Assignment** is an advanced genetic algorithm-based solution for optimally assigning workers to work centers. This project combines state-of-the-art optimization techniques with robust repair and improvement mechanisms to handle complex assignment constraints such as capacity, training levels, physical limitations, and more.

---

## Table of Contents

- [About the Project](#about-the-project)
- [Key Features](#key-features)
- [Project Structure](#project-structure)
- [Installation](#installation)
- [Usage](#usage)
- [Detailed Explanation](#detailed-explanation)
  - [1. Workcenter Database](#1-workcenter-database)
  - [2. Worker Database](#2-worker-database)
  - [3. Visualization Functions](#3-visualization-functions)
  - [4. Repair & Reassignment Mechanisms](#4-repair--reassignment-mechanisms)
  - [5. Fitness Evaluation](#5-fitness-evaluation)
  - [6. Genetic Algorithm (GA)](#6-genetic-algorithm-ga)
  - [7. Local Swap Improvement](#7-local-swap-improvement)
  - [8. Detailed Reporting & Statistics](#8-detailed-reporting--statistics)
- [Contributing](#contributing)
- [License](#license)
- [Contact](#contact)

---

## About the Project

This project is a comprehensive solution that tackles the worker assignment problem using genetic algorithms (GA). The goal is to assign workers to work centers in an optimal manner while considering multiple constraints:

- **Capacity Restrictions:** Each work center can only accommodate a fixed number of workers.
- **Training Levels:** Workers have candidate training levels for each work center, and work centers have minimum training level requirements.
- **Physical Limitations:** Some work centers require high lifting capability.
- **Additional Constraints:** Other factors such as worker limitations and forced training costs are also factored into the evaluation.

By integrating all these aspects, the project aims to deliver an optimal assignment that maximizes efficiency and minimizes costs.

---

## Key Features

- **Deterministic Random Data:**  
  A fixed random seed ensures reproducibility of experiments and consistent behavior across runs.

- **Robust Database Setup:**  
  - **Workcenter Database:** Includes production line names, capacities, training level requirements, and physical requirements.  
  - **Worker Database:** Contains worker IDs, names, limitations, high-lift capabilities, and candidate work center options with assigned training levels.

- **Advanced Repair Mechanisms:**  
  The solution incorporates a multi-step repair function that automatically fixes infeasible assignments by reassigning workers based on their candidate lists and work center capacities.

- **Dynamic Genetic Algorithm:**  
  Implements selection, two-point crossover, mutation, and even random immigrant strategies to explore the search space effectively.

- **Local Swap Improvement:**  
  After the genetic algorithm converges, a local swap improvement routine further refines the assignment by swapping workers between work centers when beneficial.

- **Detailed Reporting & Visualization:**  
  Leveraging the `rich` library, the project provides colorful and intuitive terminal outputs, including:
  - Worker and workcenter database displays.
  - Final assignment tables with detailed statuses.
  - Comprehensive statistics including overall utilization, unassigned worker counts, and forced training occurrences.
  - A dedicated report for unassigned workers with explanations and alternative suggestions.

- **Extensibility & Modularity:**  
  The code is structured in a modular way, allowing for easy modifications, testing of new genetic operators, or the integration of additional constraints.


