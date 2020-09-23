# Motor-Gearbox Database (MGDB)

The database includes this/that and is always getting bigger and we encourage users to add to it (add links), Likely want this to be its own repo hosted in the database folder. 
* Maxon
* Faulhaber 
* Tmotor 
* Allied Motion 

explanation explanation 


## Motors 
Motors are specified in files call `*_motors.csv`. When datasheets do not specify a given criteria, the value given in parenthesis should be used. A file can contain one or multiple motors but must contain the following column headers:

1. **key** - 2 letter manufacturer code followed by motor _ID. This is the *unique* identifier for each motor in the database 
1. **manufacturer** - name of company (e.g. Maxon) 
1. **ID** - product code or unique identifier e.g. 64221 or DCX10L01EBKL489)
1. **type** - "DC" or "BLDC" 
1. **V** - Nominal Voltage in Volts
1. **k_t** - Torque Constant in Nm/A
1. **R** - Winding Resistance in Ohms 
1. **L** - Inductance in Henries
1. **mass** - in kg 
1. **inertia** - kgm^2
1. **omega_nl** - No-load speed in rad/s 
1. **I_nl** - no load current in Amps (Nan)
1. **max_int_torque** - maximum intermittent torque (inf)
1. **max_int_speed** - maximum intermittent speed (inf)
1. **max_cont_speed** - maximum continuous speed (inf)
1. **max_cont_power** - maximum continuous mechanical power (inf)
1. **coulomb_friction** - in Nm (nan) 
1. **viscous_friction** - in Nms/rad (nan)
1. **Rth1** - winding-to-housing thermal resistance in K/W OR lumped thermal resistance (nan)
1. **Rth2** - housing-to-ambient thermal resistance in K/W (nan)


## Gearboxes 
Gearboxes are specified in files call `*_gearboxes.csv`. When datasheets do not specify a given criteria, the value given in parenthesis should be used. A file can contain one or multiple gearboxes but must contain the following column headers:
1. **key** - 2 letter manufacturer code followed by gearbox _ID. This is the *unique* identifier for each gearbox in the database 
1. **manufacturer** - name of company (e.g. Faulhaber)  
1. **ID** - product code or unique identifier, when a company sells multiple gear ratios with the same ID append "_R-rounded_ratio" for each unique ratio, e.g. ABCD_R83 for a ratio of 83.2:1
1. **type** - planetary, spur, cycloid, harmonic or koaxdrive
1. **stages** - number of reduction stages (1)
1. **ratio** - input to output speed ratio, a single number e.g. 90 for a 90:1 reduction. Use absolute ratios when available.  
1. **mass** - in kg (nan)
1. **inertia** - in kgm^2. This is the inertia measured at the *input* of the gearbox. 
1. **efficiency** - rated efficiency as a decimal e.g. 0.82 for 82% efficient (nan)
1. **direction** - direction, 1: output same as input, -1: direction reversed (1)
1. **max_int_torque** - maximum intermittent torque (inf) 
1. **max_cont_torque** - maximum continuous torque (inf)

## Manufacturer Codes 

* Maxon: MM 
* Faulhaber: FH
* Allied Motion: AM
* Tmotor: TM 
* Genesis: GN 

For custom designs added to the database you can use "CU" as the manufacturer code or any other 2 letter combination that is not already reserved. 

## Motor-Gearbox Compatibility 

## mgdb.m 



