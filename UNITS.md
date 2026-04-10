# Units in the GEMS Gradient Flux Pipeline

This document tracks the intended units through the oxygen and carbon dioxide gradient-flux workflow.

## Concentration

### Oxygen

- `mass_32_40`: dimensionless ratio of RGA mass 32 to argon-normalized mass 40
- `oxygen_high`, `oxygen_low`: `mL L^-1`
- `ox_high_umol_l`, `ox_low_umol_l`: `umol L^-1`

The oxygen calibration is fit in [R/calibration.R](/Users/brett/Projects/machinelab/gems/gems-processing/R/calibration.R) as:

```r
seaphox_oxygen_ml_l ~ mass_32_40
```

That means the calibrated oxygen values produced from `mass_32_40_high` and `mass_32_40_low` are in `mL L^-1`, then converted to `umol L^-1`.

### Carbon dioxide

- `mass_44_40`: dimensionless ratio of RGA mass 44 to argon-normalized mass 40
- `co2_high_umol_l`, `co2_low_umol_l`: `umol L^-1`

The CO2 calibration is fit directly in `umol L^-1`.

## Gradient

- `sensor_separation`: `m`
- `ox_gradient_umol_l_m`: `umol L^-1 m^-1`
- `co2_gradient_umol_l_m`: `umol L^-1 m^-1`

The gradient is calculated as:

```text
gradient = (high - low) / sensor_separation
```

Because `1 umol L^-1 = 1 mmol m^-3`, the gradient can also be written as:

```text
umol L^-1 m^-1 = mmol m^-4
```

So the oxygen and CO2 gradients are numerically compatible with either notation.

## Velocity

- `u`, `v`, `w`: `m s^-1`
- rotated MATLAB inputs `vx`, `vy`, `vz`: `m s^-1`

The MATLAB eddy flux script documents velocity input units as `m s^-1`.

### ADV count conversion

`gems_adv_data()` parses raw ADV `D` records and converts a subset of the integer fields to engineering units:

- `u`, `v`, and `w` velocity counts are multiplied by `0.0001`, so the output velocities are in `m/s`. The vector sends velocities as a 16 bit integer representing either 1 mm/s or 0.1mm/s, settable in configuration, and reported in the status byte (not stored here).
- `pressure` is parsed as numeric and divided by `1000`. The package treats the resulting values as engineering-unit pressure; in practice this is used as decibar-scale pressure.
- `amp1`, `amp2`, `amp3`, `corr1`, `corr2`, and `corr3` are not rescaled during parsing and remain in instrument counts / raw correlation units.

This scaling is implemented in gemstools in [R/parsers.R](R/parsers.R) in `gems_adv_data()`.


## Friction Velocity

- `Ustar`: `m s^-1`

In MATLAB, `Ustar` is derived from the covariance of turbulent velocity components:

```text
u'w' -> m^2 s^-2
Ustar = sqrt(abs(u'w')) -> m s^-1
```

## Length Scale

- `lscale`: `m`

`lscale` is a turbulent length scale parameter. In the current pipeline it is set from `length_scale = 0.4`.

## Flux

The gradient-flux equation used in R is:

```text
Flux = -Ustar * kappa * lscale * gradient
```

where:

- `Ustar`: `m s^-1`
- `kappa` (von Karman constant): dimensionless
- `lscale`: `m`
- `gradient`: `umol L^-1 m^-1`, equivalent to `mmol m^-4`

Unit expansion:

```text
(m s^-1) * (1) * (m) * (mmol m^-4) = mmol m^-2 s^-1
```

To report hourly flux:

```text
mmol m^-2 h^-1 = mmol m^-2 s^-1 * 3600
```

So the implemented hourly-rate calculation is:

```text
Flux_h = -Ustar * kappa * lscale * gradient * 3600
```

## Summary

- Concentration after calibration:
  oxygen in `mL L^-1`, then converted to `umol L^-1`
  CO2 in `umol L^-1`
- Gradient: `umol L^-1 m^-1`, equivalent to `mmol m^-4`
- Velocity: `m s^-1`
- `Ustar`: `m s^-1`
- `lscale`: `m`
- Flux before time conversion: `mmol m^-2 s^-1`
- Flux after `* 3600`: `mmol m^-2 h^-1`
