# Motor-Gearbox Database (MGDB)
This repository provides consistently formated motor/gearbox datasheets to facilitate automated motor selection for robotic applications. We hope to keep growing the list supported motors/gearboxes and welcome any contributions. As of 9/22/20 the database includes the following:
* Maxon (all motors/gearboxes)
* Faulhaber (all motors/gearboxes) 
* Tmotor (G-series Motors)
* Allied Motion (HT and Megaflux Motors)

Adding custom user defined motors and gearboxes to the database is as simple as adding a csv file as explained below. 

## Motors 
Motors are specified in files named `*_motors.csv`. When datasheets do not specify a given criteria, the value given in parenthesis should be used. A file can contain one or multiple motors but must contain the following column headers:

1. **key** - 2 letter manufacturer code (see below) followed by motor _ID. This is the *unique* identifier for each motor in the database. Keys cannot contain "\","/", or "*" 
1. **manufacturer** - name of company (e.g. Maxon) 
1. **ID** - product code or unique identifier e.g. 64221 or DCX10L01EBKL489
1. **type** - "DC" or "BLDC" 
1. **V** - Nominal Voltage in Volts
1. **k_t** - Torque Constant in Nm/A
1. **R** - Winding Resistance in Ohms 
1. **L** - Inductance in Henries
1. **mass** - in kg 
1. **inertia** - kgm<sup>2</sup>
1. **omega_nl** - No-load speed in rad/s 
1. **I_nl** - no load current in Amps (nan)
1. **I_nom** - maximum rms current in Amps (nan)
1. **max_int_torque** - maximum intermittent torque in Nm (inf)
1. **max_int_speed** - maximum intermittent speed in rad/s (inf)
1. **max_cont_speed** - maximum continuous speed in rad/s(inf)
1. **max_cont_power** - maximum continuous mechanical power in W (inf)
1. **coulomb_friction** - in Nm (nan) 
1. **viscous_friction** - in Nms/rad (nan)
1. **Rth1** - winding-to-housing thermal resistance in K/W OR lumped thermal resistance (nan)
1. **Rth2** - housing-to-ambient thermal resistance in K/W (nan)


## Gearboxes 
Gearboxes are specified in files named `*_gearboxes.csv`. When datasheets do not specify a given criteria, the value given in parenthesis should be used. A file can contain one or multiple gearboxes but must contain the following column headers:
1. **key** - 2 letter manufacturer code followed by gearbox _ID. This is the *unique* identifier for each gearbox in the database. Keys cannot contain "\","/", or "*" 
1. **manufacturer** - name of company (e.g. "Faulhaber" or "Allied_Motion")  
1. **ID** - product code or unique identifier, when a company sells multiple gear ratios with the same ID append  "_R rounded_ratio" for each unique ratio, e.g. ABCD_R83 for a ratio of 83.2:1
1. **type** - planetary, spur, cycloid, harmonic or koaxdrive
1. **stages** - number of reduction stages (1)
1. **ratio** - input to output speed ratio, a single number e.g. 90 for a 90:1 reduction. Use absolute ratios when available. 
1. **mass** - in kg (nan)
1. **inertia** - in kgm<sup>2</sup>. This is the inertia measured at the *input* of the gearbox. 
1. **efficiency** - rated efficiency as a decimal e.g. 0.82 for 82% efficient (nan)
1. **direction** - direction, 1: output same as input, -1: direction reversed (1)
1. **max_int_torque** - maximum intermittent torque in Nm (inf) 
1. **max_cont_torque** - maximum continuous torque in Nm (inf)

## Manufacturer Codes 

* Maxon: MM 
* Faulhaber: FH
* Allied Motion: AM
* Tmotor: TM 
* Genesis: GN 

For custom designs added to the database you can use "CU" as the manufacturer code or any other 2 letter combination that is not already reserved. 

## Motor-Gearbox Compatibility 
Compatibility between motors and gearboxes is added using files named  `*_compatibility.csv`. The first entry of each row of the file gives the motor key and all subsequent entries in that row give the keys of compatible gearboxes. Gearbox keys in compatibility files can be *partial keys*, for example "FH_26A_R*" is a partial key that denotes *all* gearboxes whose keys start with "FH_26_A_R". 

## mgdb.m 
`mdgb.m` provides a Matlab class for buidling a database from motor/gearbox/compatibility files. Type `doc mgdb` at the Matlab command prompt for documentation.


