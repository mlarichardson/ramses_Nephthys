Star-Formation and Feedback Patch
=================================

Star-Formation
--------------

### Overview ###

A new star-formation criterion has been implemented from the patch sent by
Julien Devriendt and Taysun Kimm. It uses a "virial parameter" criterion, to
form stars only where the gas is Jeans unstable (taking into account the
small-scale turbulence via a "turbulent Jeans length"). The reference paper
for this criterion is Federrath & Klessen 2012, ApJ 761:156.

For the star-forming regions, the star-formation timescale is redefined in the
star_formation.f90 routine, once again following the formalism of F&K 2012.

This criterion is triggered by the sf_hopkins flag in the PHYSICS namelist.

### Main changes ###

*   __star_formation.f90__: main implementation of the SF criterion. The selection
    of new star-forming cells is changed, as well as the PoissMean and tstar:

       PoissMean = dtnew(ilevel)*sfr_ff(i)/tstar*mcell/mstar

*   __read_hydro_params.f90__: added "sf_hopkins" in the namelist

*   __amr_parameter.f90__: added "sf_hopkins" logical


Feedback
--------

### Overview ###

A new SN feedback recipe has been implemented by Taysun Kimm, as described in
the appendix A of Kimm & Cen 2014, ApJ 788:121. The idea is to inject the
energy either in a momentum or energy conserving way, depending on which stage
of the Sedov expansion phase is supposed to be described by the simulation.

### Main changes ###

*   __feedback.f90__: a new kinetic_feedback_cell has been added, along with some
    convenience routines. There is no free parameter (except the delay before
    the stars go SN).

*   __amr_step.f90__: In the "Kinetic feedback from giant molecular clouds" part,
    call to kinetic_feedback_cell instead of kinetic_feedback


Minor changes
-------------

*   Following Taysun Kimm's version, a small change has been made to
    __newdt_fine.f90__ to get the outputs exactly when we want.

*   ~~Slight change in __rho_fine.f90__ for something sink related that seemed
    to be not working (maybe due to the work in progress on the sinks in the
    main ramses version?)~~
