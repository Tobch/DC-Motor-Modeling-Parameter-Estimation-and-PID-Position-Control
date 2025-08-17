# DC-Motor-Modeling-Parameter-Estimation-and-PID-Position-Control


# DC Motor Modeling and Control -- README

## 📌 Project Overview

This project models, identifies, and controls a DC motor using
MATLAB/Simulink. The workflow includes system modeling, parameter
estimation from measured data, transfer function (TF) identification,
and PID controller design.

------------------------------------------------------------------------

## ⚙️ Steps Performed

### 1. System Modeling

-   The DC motor was represented as a **second-order system**.\
-   Transfer function form:\
    \[ G(s) =
    `\frac{b_2 s^2 + b_1 s + b_0}{a_2 s^2 + a_1 s + a_0}`{=tex} \]

### 2. Parameter Estimation

-   Experimental data was logged using `read_dc_motor` and
    `log_motor_sequence`.\
-   MATLAB's **Parameter Estimation app** identified updated transfer
    functions from measured input--output data.\
-   Example identified TF (speed response):\
    \[ G(s) =
    `\frac{-1642.8s^2 + 2.87 \times 10^5 s + 1.39 \times 10^6}{332.3s^2 + 5824.2s + 17198.3}`{=tex}
    \]

### 3. Control Analysis

-   **Poles** were extracted from the denominator.\
-   **Damping ratio** calculated as:\
    \[ `\zeta `{=tex}= `\frac{a_1/a_2}{2\sqrt{a_0/a_2}}`{=tex} \]
-   Example results:
    -   ωₙ = 7.19 rad/s\
    -   ζ = 1.22 → **Overdamped** system\
    -   Poles: -13.77, -3.76 (real, negative)

### 4. PID Controller Design

-   PID controller (`Kp`, `Ki`, `Kd`) designed using MATLAB (`pidtune`)
    with the **updated TF**.\
-   Controller implemented in **Simulink** for motor position control.

------------------------------------------------------------------------

## 📊 Results

-   Simulated responses closely matched measured motor responses after
    parameter estimation.\
-   Overdamped behavior confirmed (no oscillations, slow response).\
-   PID controller improved rise time and settling time while
    maintaining stability.

------------------------------------------------------------------------

## 📁 Deliverables

-   **Simulink models**: `dc_motor_sim.slx`, `motor_pid_control.slx`\
-   **Report** (DOCX) -- includes methodology, results, and discussion.\
-   **This README** -- summary of workflow.

------------------------------------------------------------------------

