"""
@authors:   Luis Juarez Mercedes
            Pedro Vidal Arias
"""

"""
Enable Actuators
"""

from adafruit_servokit import ServoKit

motor = ServoKit(channels=16)
# Set pulse width range
motor.servo[0].set_pulse_width_range(600, 2400) # Gripper
motor.servo[1].set_pulse_width_range(500, 2500) # Axis1
motor.servo[2].set_pulse_width_range(500, 2500) # Axis2
motor.servo[3].set_pulse_width_range(500, 2500) # Axis3
motor.servo[4].set_pulse_width_range(500, 2500) # Axis4
motor.servo[5].set_pulse_width_range(500, 2500) # Axis5
motor.servo[6].set_pulse_width_range(500, 2500) # Axis6
# Set actuation range
motor.servo[0].actuation_range = 180    # Gripper
motor.servo[1].actuation_range = 270    # Axis1
motor.servo[2].actuation_range = 270    # Axis2
motor.servo[3].actuation_range = 270    # Axis3
motor.servo[4].actuation_range = 270    # Axis4
motor.servo[5].actuation_range = 270    # Axis5
motor.servo[6].actuation_range = 270    # Axis6


