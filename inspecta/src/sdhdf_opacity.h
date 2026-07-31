//  Copyright (C) 2019, 2020, 2021, 2022 George Hobbs

/*
 *    This file is part of INSPECTA.
 *
 *    INSPECTA is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *    INSPECTA is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *    You should have received a copy of the GNU General Public License
 *    along with INSPECTA.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef __SDHDF_OPACITY_H
#define __SDHDF_OPACITY_H

// ITU-R P.676 gaseous attenuation -- the table-free approximate (Annex 2 Sec 1)
// slant-path method, driven by ground-based weather (temperature, pressure,
// relative humidity). Ported from atnf-parkes-opacity's P676.java, itself a
// direct port of the itur Python package's P.676-10 approximate method, so
// values match that reference. Valid for f <= 350 GHz and elevation 5-90 deg.

// Water-vapour density (g/m^3) from temperature (deg C) and relative humidity
// (%), via the Arden Buck saturation-pressure equation, clipped to [0.1,30].
double sdhdf_opacity_waterVapourDensity(double tempC,double relHumidityPct);

// Slant-path gaseous attenuation (dB) at elevation (deg, 90 = zenith).
double sdhdf_opacity_attenuationDb(double freqGHz,double elevationDeg,double rho,double pressureHpa,double tempK);

// Slant-path opacity tau (nepers) at elevation (deg).
double sdhdf_opacity_tau(double freqGHz,double elevationDeg,double rho,double pressureHpa,double tempK);

// Multiplicative correction factor to apply to calibrated flux/antenna
// temperature to compensate for atmospheric gaseous absorption:
// exp(tau_slant), clamped to [1.0,3.0] as per the reference implementation
// (atnf-pkspoint pyspot.py). freqMHz is the observing frequency; tempC,
// pressureHpa, relHumidityPct are ground weather at the antenna.
double sdhdf_opacity_correctionFactor(double freqMHz,double elevationDeg,double tempC,double pressureHpa,double relHumidityPct);

#endif
