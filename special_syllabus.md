# Special Syllabus Agreement (Spesialpensum) – NMBU

---

## 1. Contract Between Student and Academic Responsible

```text
Main supervisor for the special syllabus: Jon Glenn Gjevestad
Co-supervisor for the special syllabus: Narve Kjørsvik
```
---

## 2. Information About the Special Syllabus

### Norwegian Title
```text
Anvendt treghetsnavigasjon
```

### English Title
```text
Applied inertial navigation
```

---

### Learning Goals
```text
The overarching goal of the special syllabus is to provide the student with a solid theoretical understanding and practical competence in strapdown inertial navigation systems (INS) aided by satellite navigation (GNSS). 

The student will implement a functional processing framework in GNU Octave (or MATLAB) for processing raw measurements from an Inertial Measurement Unit (IMU) loosely coupled with GNSS positions via an error-state (complementary) Kalman filter. Using this framework, the student will process real-world high-dynamics platform data (e.g., aircraft or helicopter), evaluate system performance, and compare estimated trajectories against reference solutions from commercial post-processing software.
```

---

### Learning Outcomes
```text
Upon successful completion of the special syllabus, the candidate will have achieved the following learning outcomes:

Knowledge:
• Explain the fundamental operating principles and sensor mechanics of strapdown inertial measurement units (IMUs), including gyroscopes and accelerometers.
• Formulate the kinematic equations of motion across relevant geodetic and navigation coordinate reference frames (ECEF, NED, body frame) and understand attitude representations (Euler angles, direction cosine matrices, quaternions).
• Understand error sources and stochastic error modeling for inertial sensors (biases, scale factors, drift, random walk).
• Describe the theory and architecture of complementary / error-state extended Kalman filtering for loose GNSS/INS integration.

Skills:
• Implement a modular computational framework in GNU Octave to read and process high-rate IMU accelerometer and gyroscope observations.
• Develop and tune a discrete-time complementary Kalman filter to estimate position, velocity, attitude errors, and sensor biases using aiding GNSS position fixes.
• Process and analyze a real high-dynamics kinematic flight dataset with appropriate initial alignment and trajectory propagation.
• Validate filter outputs by conducting consistency checks, innovation analysis, and benchmark comparisons against commercial post-processing navigation software.

General Competence:
• Critically assess the accuracy and limitations of GNSS aided multisensor navigation systems.
```

---

### Credits, Grading, and Assessment Form
```text
Number of credits: 5
Select the marking system for your special syllabus: Pass/Fail (Bestått / Ikke bestått)
Form of assessment: Report
```

---

#### Description of Assessment Form / Deliverable
```text
The special syllabus is assessed on a Pass/Fail (Bestått / Ikke bestått) basis evaluated on an individual comprehensive technical report. The report must thoroughly document:
1. Mathematical foundation and system architecture (mechanization equations, error-state model, and filter equations).
2. Description of the GNU Octave program code, algorithmic structure, and modular implementation.
3. Analysis and discussion of empirical results from processing the real high-dynamics test dataset, including trajectory plots, attitude states, covariance estimates, and comparative evaluation against commercial benchmark software.
```

---

### Period and Submission Deadlines
```text
Completion of the special syllabus (Year and block/parallel): 2026 autumn parallel
Submission deadline: 15.12.2026
Name of proposed external sensor: Kjetil Bergh Aanonsen (FFI)
```

---

## 3. Attachment: Framework, Content, Progress Plan & Literature Reference List

### 3.1 Course Content and Scope (5 ECTS)
A 5 ECTS special syllabus corresponds to approximately 125–150 hours of dedicated student workload. The syllabus centers on strapdown inertial navigation systems, multisensor integration, and Kalman filtering:
1. **Coordinate Systems & Kinematics**: Inertial frame, Earth-Centered Earth-Fixed (ECEF), North-East-Down (NED), and Body frame. Attitude representations and transformations.
2. **Inertial Sensors & Mechanization**: Gyroscope and accelerometer observation equations, Earth rate and gravity compensation, strapdown integration algorithms.
3. **Stochastic Modeling & Optimal Estimation**: Continuous and discrete state-space models, white noise, random walk, Gauss-Markov processes, and complementary / error-state Extended Kalman Filter (EKF).
4. **Loosely Coupled GNSS/INS Integration**: Measurement update using GNSS positions, correction of position, velocity, attitude, and sensor biases.
5. **Practical Implementation**: Scripting the processing pipeline in GNU Octave, processing high-dynamics airborne data, validating against benchmark commercial software.

---

### 3.2 Literature Reference List (Pensumliste)

**Primary Textbook:**
* **Farrell, Jay A. (2008).** *Aided Navigation: GPS with High Rate Sensors*. McGraw-Hill Education. ISBN: 9780071493291.
  * *Selected Chapters (approx. 200 pages):*
    * **Chapter 1:** Overview (Introduction to aided navigation systems)
    * **Chapter 2:** Reference Frames & Coordinate Systems (Coordinate frames, transformations, attitude representations)
    * **Chapter 3:** Deterministic Systems (State-space representations, kinematics)
    * **Chapter 4:** Stochastic Processes (Random variables, correlation, power spectral density, sensor noise models)
    * **Chapter 5:** Optimal State Estimation (Kalman filter formulation, discrete-time propagation and measurement updates)
    * **Chapter 7:** Navigation System Design (Error-state formulation, complementary filtering architecture)
    * **Chapter 11:** Aided Inertial Navigation (Strapdown mechanization, loosely coupled GNSS/INS integration, bias tracking)

**Supplementary Materials:**
* Compendium and lecture notes on strapdown inertial navigation and Kalman filtering provided by the academic supervisors (NMBU).
* Documentation and interface specifications for the provided high-dynamics flight dataset and commercial benchmark software.

---
