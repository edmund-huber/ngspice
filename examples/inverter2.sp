Inverter example circuit
* This netlist demonstrates the following:
* global nodes (vdd, gnd)
* autostop (.tran defines simulation end as 4ns but simulation stops at
* 142.5ps when .measure statements are evaluated)
* scale (all device units are in microns)
* model binning (look in device.values file for which bin chosen)
*
*     m.x1.mn:
*     model              = nch.2
*
*     m.x1.mp:
*     model              = pch.2
*
* parameters
* parameterized subckt
* vsrc with repeat
* .measure statements for delay and an example ternary operator
* Additional .measure statements to the previous example
* to measure the min, avg, rms and integ values for the VDD voltage
* source.  These .measure statements disable 'autostop'.
* device listing and parameter listing
* You can run the example circuit with this command:
*
* ngspice inverter2.sp


* global nodes
.global vdd gnd

* autostop -- stop simulation early if .measure statements done
* scale    -- define scale factor for mosfet device parameters (l,w,area,perimeter)
.option autostop
.option scale = 1e-6

* model binning
.model nch.1 nmos ( version=4.4 level=54 lmin=0.1u lmax=20u wmin=0.1u wmax=10u  )
.model nch.2 nmos ( version=4.4 level=54 lmin=0.1u lmax=20u wmin=10u  wmax=100u )
.model pch.1 pmos ( version=4.4 level=54 lmin=0.1u lmax=20u wmin=0.1u wmax=10u  )
.model pch.2 pmos ( version=4.4 level=54 lmin=0.1u lmax=20u wmin=10u  wmax=100u )

* parameters
.param vp     = 1.0v
.param lmin   = 0.10
.param wmin   = 0.12
.param plmin  = 'lmin'
.param nlmin  = 'lmin'
.param wpmin  = 'wmin'
.param wnmin  = 'wmin'
.param drise  = 400ps
.param dfall  = 100ps
.param trise  = 100ps
.param tfall  = 100ps
.param period = 1ns

* parameterized subckt
.subckt inv in out pw='wpmin' pl='plmin' nw='wnmin' nl='nlmin'
mp out in vdd vdd pch w='pw' l='pl'
mn out in gnd gnd nch w='nw' l='nl'
.ends

v0 vdd gnd 'vp'

* vsrc with repeat
v1 in gnd pwl
+ 0ns                       'vp'
+ 'dfall-0.8*tfall'         'vp'
+ 'dfall-0.4*tfall'         '0.9*vp'
+ 'dfall+0.4*tfall'         '0.1*vp'
+ 'dfall+0.8*tfall'         0v
+ 'drise-0.8*trise'         0v
+ 'drise-0.4*trise'         '0.1*vp'
+ 'drise+0.4*trise'         '0.9*vp'
+ 'drise+0.8*trise'         'vp'
+ 'period+dfall-0.8*tfall'  'vp'
+ r='dfall-0.8*tfall'

x1 in out inv pw=60 nw=20
c1 out gnd 220fF

.tran 1ps 4ns

.measure tran inv_delay trig v(in)  val='vp/2'   fall=1 targ v(out) val='vp/2'   rise=1
.measure tran out_slew  trig v(out) val='0.2*vp' rise=1 targ v(out) val='0.8*vp' rise=1
.measure tran delay_chk param='(inv_delay < 100ps) ? 1 : 0'
.measure tran v0_min    min   i(v0) from='dfall' to='dfall+period'
.measure tran v0_avg    avg   i(v0) from='dfall' to='dfall+period'
.measure tran v0_integ  integ i(v0) from='dfall' to='dfall+period'
.measure tran v0_rms    rms   i(v0) from='dfall' to='dfall+period'

.control
run
show    *     > device.values
listing param > param.values
.endc

.end
