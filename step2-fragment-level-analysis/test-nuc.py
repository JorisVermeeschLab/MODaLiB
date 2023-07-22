import time
from typing import Tuple

import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from modalib.core import Sample


class PIDController:

    def __init__(self):
        self.kp = 0.005
        self.ki = 0.0001
        self.kd = 0.005
        self.previous_error = 0.0
        self.integral = 0.0

    def update(self, setpoint, measured_value):
        error = setpoint - measured_value
        proportional = error
        dt = 1  # TODO
        self.integral += error * dt
        derivative = (error - self.previous_error) / dt
        output = self.kp * proportional + self.ki * self.integral + self.kd * derivative
        self.previous_error = error
        return output


class Oscillator:

    PHASE_INCREMENT = 0.015
    PERIOD_INCREMENT = 0.015
    MIN_PERIOD = 120.0
    MAX_PERIOD = 200.0

    def __init__(self):
        self.phi0 = 0.0
        self.period = 167.0
        self.pid = PIDController()

    def __call__(self, genomic_coord):
        return np.cos(self.angle(genomic_coord))

    def derivative(self, genomic_coord):
        return -np.sin(self.angle(genomic_coord))

    def angle(self, genomic_coord):
        return 2. * np.pi * genomic_coord / self.period + self.phi0

    def increase_period(self):
        self.period += Oscillator.PERIOD_INCREMENT
        self.period = np.clip(self.period, Oscillator.MIN_PERIOD, Oscillator.MAX_PERIOD)

    def decrease_period(self):
        self.period -= Oscillator.PERIOD_INCREMENT
        self.period = np.clip(self.period, Oscillator.MIN_PERIOD, Oscillator.MAX_PERIOD)

    def increase_phase(self):
        self.phi0 += Oscillator.PHASE_INCREMENT

    def decrease_phase(self):
        self.phi0 -= Oscillator.PHASE_INCREMENT

    def update_new(self, start, end):
        y1 = self.__call__(start)
        derivative1 = self.derivative(start)
        y2 = self.__call__(end)
        derivative2 = self.derivative(end)
        if (y1 > 0) and (y2 > 0):
            if (derivative1 > 0) and (derivative2 > 0):
                self.decrease_phase()
            elif (derivative1 < 0) and (derivative2 < 0):
                self.increase_phase()
            if end - start > self.period:
                self.increase_period()
            else:
                self.decrease_period()

    def update(self, start, end):
        y1 = self.__call__(start)
        derivative1 = self.derivative(start)
        y2 = self.__call__(end)
        derivative2 = self.derivative(end)
        if (y1 > 0) and (y2 > 0):
            if (derivative1 > 0) and (derivative2 > 0):
                self.decrease_phase()
            elif (derivative1 < 0) and (derivative2 < 0):
                self.increase_phase()
            if end - start > self.period:
                self.increase_period()
            else:
                self.decrease_period()


sample_names = [
    'GC131877', 'GC116736', 'GC131879', 'GC116728', 'GC110691'
]

samples = [Sample(f'../step1-methyldackel/out/{sample_name}') for sample_name in sample_names]


def trajectory(sample: Sample) -> Tuple[np.ndarray, np.ndarray]:
    osc = Oscillator()
    xs, ys = [], []
    count = 0
    score = 0.
    for t, read in enumerate(sample):
        if read.chr_id > 1:
            break
        start = read.start
        end = read.start + read.size
        score += -osc(start)
        score += -osc(end)
        count += 2

        osc.update(start, end)
        xs.append(start)
        ys.append(osc.period)

    score /= count

    print(f'Score: {score}')
    return np.asarray(xs), np.asarray(ys)


xs, ys = trajectory(samples[0])
plt.plot(xs, ys)
xs, ys = trajectory(samples[1])
plt.plot(xs, ys)
xs, ys = trajectory(samples[2])
plt.plot(xs, ys)
xs, ys = trajectory(samples[3])
plt.plot(xs, ys)
xs, ys = trajectory(samples[4])
plt.plot(xs, ys)
plt.savefig('traj.png')

print('Finished')
