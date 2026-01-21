"""
Eddy Flux Batch Processing Script

Calculates eddy covariance fluxes from acoustic Doppler velocity and biogeochemical data.

Input:
    - data file with columns: time (h), vx, vy, vz (m/s), o2 (umol/L), 
      pH (units), pressure (dbar), SNR, correlation
    - timestamp file with burst begin/end times

Output:
    - Burst mean fluxes and parameters
    - Cumulative fluxes (optional)

Date created: 6/5/2013; updated 9/14/18; modified for batch input 11/21/24
Translated to Python: 12/2025

UNTESTED!
Not currenty used for flux processing
"""

import numpy as np
import pandas as pd
from scipy import signal
from pathlib import Path
import sys


class EddyFluxProcessor:
    """Process eddy covariance data and calculate fluxes."""
    
    def __init__(self, infile, timestampfile, flagrotate=1, rotation_index=None):
        """
        Initialize the eddy flux processor.
        
        Parameters
        ----------
        infile : str
            Input data filename with velocity, O2, pH, pressure data
        timestampfile : str
            File with burst begin/end timestamps
        flagrotate : int
            0=no rotation, 1=standard rotation, 2=planar rotation (default=1)
        rotation_index : int
            Index for planar rotation angles if flagrotate==2
        """
        self.infile = infile
        self.timestampfile = timestampfile
        self.flagrotate = flagrotate
        
        # Output filenames
        self.outfile2 = f"outfile2_{Path(infile).name}"
        self.outfile3 = f"outfile3_{Path(infile).name}"
        
        # Processing parameters
        self.mheight = 0.6  # Measuring height (m)
        self.heading = 0.0  # Vector heading (degrees)
        self.hz = 4  # Sampling frequency (Hz)
        self.windowsize = self.hz * 60 * 5 + 1
        self.window = self.hz * 60 * 14 + 1  # Hamming window size
        self.Td = 5  # Wave frequency cutoff period (s)
        self.B = 4  # Bursts per hour
        self.dist = 0.025  # Distance between sensors (m)
        self.flowlag = 0.5  # Sensor response time lag (s)
        self.flagshift = 1  # Apply timeshift correction
        self.flagwv = 0.045  # Min velocity threshold for Cd/phi2 (m/s)
        self.flagwrite = 0  # Write cumulative fluxes (slow if True)
        
        # Planar rotation angles
        self.planar_rotations = np.array([
            0.00, 0.00, 0.00, -2.65, -2.59, 0.70, 2.31, -0.20, 
            1.16, 0.81, 0.90, -3.95, 1.49, 0.58, -4.44
        ])
        if rotation_index is not None and flagrotate == 2:
            self.planarzrot = np.radians(self.planar_rotations[rotation_index])
        else:
            self.planarzrot = np.radians(-5.5725)
    
    def load_data(self):
        """Load input data and timestamps."""
        print("Importing data...")
        
        # Load data file (skip header)
        data = pd.read_csv(self.infile, sep=r'\s+', skiprows=1, header=None, dtype=float)
        
        self.time = data.iloc[:, 0].values.astype(float)
        self.vx = data.iloc[:, 1].values.astype(float)
        self.vy = data.iloc[:, 2].values.astype(float)
        self.vz = data.iloc[:, 3].values.astype(float)
        self.o2 = data.iloc[:, 4].values.astype(float)
        self.ph = data.iloc[:, 5].values.astype(float)
        self.pres = data.iloc[:, 6].values.astype(float)
        self.SNR = data.iloc[:, 7].values.astype(float)
        self.correlation = data.iloc[:, 8].values.astype(float)
        
        # Load timestamps
        self.timestamps = pd.read_csv(self.timestampfile, sep=r'\s+', header=None, dtype=float).values
        self.n_bursts = len(self.timestamps)
    
    def find_burst_indices(self):
        """Find time indices for burst start/end."""
        self.j1 = np.zeros(self.n_bursts, dtype=int)
        self.j2 = np.zeros(self.n_bursts, dtype=int)
        
        for i in range(self.n_bursts):
            self.j1[i] = np.where(self.time > self.timestamps[i, 0])[0][0]
            self.j2[i] = np.where(self.time < self.timestamps[i, 1])[0][-1]
    
    def convert_ph(self):
        """Convert pH to hydrogen ion concentration."""
        self.H = (10.0 ** (-self.ph)) * 1e6  # uMol L-1
        
        # Calculate storage terms (6-hour moving average)
        storage_window = int(self.hz * 3600 * 6)
        self.o2_storage = np.full_like(self.o2, np.nan)
        self.ph_storage = np.full_like(self.H, np.nan)
        
        for i in range(self.n_bursts):
            idx_start = self.j1[0]
            idx_end = self.j2[-1]
            self.o2_storage[idx_start:idx_end] = self._moving_mean(
                self.o2[idx_start:idx_end], storage_window
            )
            self.ph_storage[idx_start:idx_end] = self._moving_mean(
                self.H[idx_start:idx_end], storage_window
            )
    
    @staticmethod
    def _moving_mean(data, window):
        """Calculate moving mean."""
        return pd.Series(data).rolling(window=window, center=True).mean().values
    
    def calculate_burst_means(self):
        """Calculate mean values for each burst."""
        print("Calculating burst means...")
        
        self.timemean = np.zeros(self.n_bursts)
        self.vxmean = np.zeros(self.n_bursts)
        self.vymean = np.zeros(self.n_bursts)
        self.vzmean = np.zeros(self.n_bursts)
        self.vmean = np.zeros(self.n_bursts)
        self.o2mean = np.zeros(self.n_bursts)
        self.phmean = np.zeros(self.n_bursts)
        self.Hmean = np.zeros(self.n_bursts)
        self.presmean = np.zeros(self.n_bursts)
        self.Hs = np.zeros(self.n_bursts)  # Wave height
        self.SNRmean = np.zeros(self.n_bursts)
        self.correlationmean = np.zeros(self.n_bursts)
        
        for i in range(self.n_bursts):
            idx = slice(self.j1[i], self.j2[i] + 1)
            
            self.timemean[i] = np.mean(self.time[idx])
            self.vxmean[i] = np.mean(self.vx[idx])
            self.vymean[i] = np.mean(self.vy[idx])
            self.vzmean[i] = np.mean(self.vz[idx])
            self.vmean[i] = np.sqrt(self.vxmean[i]**2 + self.vymean[i]**2 + self.vzmean[i]**2)
            self.o2mean[i] = np.mean(self.o2[idx])
            self.phmean[i] = np.mean(self.ph[idx])
            self.Hmean[i] = np.mean(self.H[idx])
            self.presmean[i] = np.mean(self.pres[idx])
            self.Hs[i] = 4.0 * np.std(self.pres[idx])
            self.SNRmean[i] = np.mean(self.SNR[idx])
            self.correlationmean[i] = np.mean(self.correlation[idx])
    
    def calculate_rotations(self):
        """Calculate coordinate rotations."""
        print("Calculating rotations...")
        
        n_data = len(self.time)
        self.theta = np.zeros(self.n_bursts)
        self.phi = np.zeros(self.n_bursts)
        self.vrot2 = np.zeros((n_data, 3))
        self.vxrot = np.zeros(self.n_bursts)
        self.vyrot = np.zeros(self.n_bursts)
        self.vzrot = np.zeros(self.n_bursts)
        self.thetad = np.zeros(self.n_bursts)
        self.phid = np.zeros(self.n_bursts)
        
        for i in range(self.n_bursts):
            idx = slice(self.j1[i], self.j2[i] + 1)
            
            # First rotation: align horizontal velocity with x-axis
            self.theta[i] = np.arctan2(self.vymean[i], self.vxmean[i])
            rot1 = self._rotation_matrix_z(self.theta[i])
            burst_data = np.column_stack([self.vx[idx], self.vy[idx], self.vz[idx]])
            vrot1 = burst_data @ rot1
            
            # Second rotation: align vertical velocity with z-axis
            if self.flagrotate == 1:
                self.phi[i] = np.arctan2(vrot1[:, 2].mean(), vrot1[:, 0].mean())
            elif self.flagrotate == 2:
                self.phi[i] = self.planarzrot
            
            rot2 = self._rotation_matrix_y(self.phi[i])
            vrot2_burst = vrot1 @ rot2
            self.vrot2[self.j1[i]:self.j2[i]+1] = vrot2_burst
            
            self.vxrot[i] = vrot2_burst[:, 0].mean()
            self.vyrot[i] = vrot2_burst[:, 1].mean()
            self.vzrot[i] = vrot2_burst[:, 2].mean()
            
            self.thetad[i] = np.degrees(self.theta[i])
            self.phid[i] = np.degrees(self.phi[i])
    
    @staticmethod
    def _rotation_matrix_z(angle):
        """Z-axis rotation matrix."""
        cos_a = np.cos(angle)
        sin_a = np.sin(angle)
        return np.array([
            [cos_a, -sin_a, 0],
            [sin_a, cos_a, 0],
            [0, 0, 1]
        ])
    
    @staticmethod
    def _rotation_matrix_y(angle):
        """Y-axis rotation matrix."""
        cos_a = np.cos(angle)
        sin_a = np.sin(angle)
        return np.array([
            [cos_a, 0, -sin_a],
            [0, 1, 0],
            [sin_a, 0, cos_a]
        ])
    
    def calculate_horizontal_angle(self):
        """Calculate horizontal current direction."""
        self.horz = np.zeros(self.n_bursts)
        
        for i in range(self.n_bursts):
            thetad2 = self.thetad[i] + self.heading
            
            # Normalize to 0-360 range
            while thetad2 < 0:
                thetad2 += 360
            while thetad2 > 360:
                thetad2 -= 360
            
            self.horz[i] = thetad2
    
    def apply_timeshift(self):
        """Apply timeshift correction based on flow lag and sensor spacing."""
        self.shift = np.zeros(self.n_bursts, dtype=int)
        
        for i in range(self.n_bursts):
            self.shift[i] = int(np.floor(self.dist / self.vmean[i] * self.hz + 
                                         self.flowlag * self.hz))
        
        if self.flagshift == 1:
            o2_old = self.o2.copy()
            H_old = self.H.copy()
            self.o2 = np.full_like(self.o2, np.nan)
            self.H = np.full_like(self.H, np.nan)
            
            for i in range(self.n_bursts):
                burst_start = self.j1[i]
                burst_end = self.j2[i] + 1
                src_start = self.j1[i] + self.shift[i]
                src_end = min(self.j2[i] + self.shift[i] + 1, len(o2_old))
                
                # Only copy if within bounds
                dest_len = burst_end - burst_start
                src_len = src_end - src_start
                copy_len = min(dest_len, src_len)
                
                if copy_len > 0:
                    self.o2[burst_start:burst_start+copy_len] = o2_old[src_start:src_start+copy_len]
                    self.H[burst_start:burst_start+copy_len] = H_old[src_start:src_start+copy_len]
    
    def calculate_flux2(self):
        """Calculate flux2 using linear detrending."""
        print("Calculating flux2 (linear detrending)...")
        
        self.flux2o2 = np.zeros(self.n_bursts)
        self.flux2ph = np.zeros(self.n_bursts)
        self.TKE2 = np.zeros(self.n_bursts)
        self.wv2 = np.zeros(self.n_bursts)
        
        for i in range(self.n_bursts):
            idx = slice(self.j1[i], self.j2[i] + 1)
            
            # Choose velocity data based on rotation
            if self.flagrotate == 0:
                vx_burst = self.vx[idx]
                vy_burst = self.vy[idx]
                vz_burst = self.vz[idx]
            else:
                vx_burst = self.vrot2[idx, 0]
                vy_burst = self.vrot2[idx, 1]
                vz_burst = self.vrot2[idx, 2]
            
            # Linear fit and residuals
            t_burst = self.time[idx]
            o2fit = np.polyfit(t_burst, self.o2[idx], 1)
            phfit = np.polyfit(t_burst, self.H[idx], 1)
            vxfit = np.polyfit(t_burst, vx_burst, 1)
            vyfit = np.polyfit(t_burst, vy_burst, 1)
            vzfit = np.polyfit(t_burst, vz_burst, 1)
            
            o2_prime = self.o2[idx] - np.polyval(o2fit, t_burst)
            ph_prime = self.H[idx] - np.polyval(phfit, t_burst)
            vx_prime = vx_burst - np.polyval(vxfit, t_burst)
            vy_prime = vy_burst - np.polyval(vyfit, t_burst)
            vz_prime = vz_burst - np.polyval(vzfit, t_burst)
            
            # Fluxes (convert from h to h-1)
            self.flux2o2[i] = np.mean(vz_prime * o2_prime) * 3600
            self.flux2ph[i] = np.mean(vz_prime * ph_prime) * 3600
            
            # Turbulent kinetic energy
            self.TKE2[i] = np.sqrt(np.mean(vx_prime**2) + np.mean(vy_prime**2) + 
                                   np.mean(vz_prime**2))
            
            # Wave velocity (horizontal)
            self.wv2[i] = np.sqrt(np.mean((self.vx[idx] - self.vxmean[i])**2) + 
                                 np.mean((self.vy[idx] - self.vymean[i])**2))
    
    def calculate_flux3(self):
        """Calculate flux3 using high-pass filtered data."""
        print("Calculating flux3 (high-pass filtering)...")
        
        self.flux3o2 = np.zeros(self.n_bursts)
        self.flux3ph = np.zeros(self.n_bursts)
        self.flux3o2_stor = np.zeros(self.n_bursts)
        self.flux3ph_stor = np.zeros(self.n_bursts)
        self.TKE3 = np.zeros(self.n_bursts)
        self.wv3 = np.zeros(self.n_bursts)
        self.Flux3cpsd = np.zeros(self.n_bursts)
        self.Flux3cpsd_low = np.zeros(self.n_bursts)
        self.Flux3cpsd_ph = np.zeros(self.n_bursts)
        self.Flux3cpsd_ph_low = np.zeros(self.n_bursts)
        self.Flux3cpsd_low_stor = np.zeros(self.n_bursts)
        self.Flux3cpsd_ph_low_stor = np.zeros(self.n_bursts)
        self.Flux_cpsd_xy = np.zeros(self.n_bursts)
        self.Flux_cpsd_xy_low = np.zeros(self.n_bursts)
        self.Ustar = np.zeros(self.n_bursts)
        self.Zzero = np.zeros(self.n_bursts)
        self.EddyDiff = np.zeros(self.n_bursts)
        self.Cd = np.full(self.n_bursts, np.nan)
        self.phi2 = np.full(self.n_bursts, np.nan)
        
        for i in range(self.n_bursts):
            idx = slice(self.j1[i], self.j2[i] + 1)
            
            # Choose velocity data based on rotation
            if self.flagrotate == 0:
                vx_burst = self.vx[idx]
                vy_burst = self.vy[idx]
                vz_burst = self.vz[idx]
            else:
                vx_burst = self.vrot2[idx, 0]
                vy_burst = self.vrot2[idx, 1]
                vz_burst = self.vrot2[idx, 2]
            
            # High-pass filter using moving mean subtraction
            vx_prime = vx_burst - self._moving_mean(vx_burst, self.windowsize)
            vy_prime = vy_burst - self._moving_mean(vy_burst, self.windowsize)
            vz_prime = vz_burst - self._moving_mean(vz_burst, self.windowsize)
            o2_prime = self.o2[idx] - self._moving_mean(self.o2[idx], self.windowsize)
            ph_prime = self.H[idx] - self._moving_mean(self.H[idx], self.windowsize)
            
            # Fluxes
            self.flux3o2[i] = np.mean(vz_prime * o2_prime) * 3600
            self.flux3ph[i] = np.mean(vz_prime * ph_prime) * 3600
            
            # Storage flux
            o2_start = np.mean(self.o2_storage[self.j1[i]:self.j1[i]+self.hz*60])
            o2_end = np.mean(self.o2_storage[self.j2[i]-self.hz*60:self.j2[i]])
            ph_start = np.mean(self.ph_storage[self.j1[i]:self.j1[i]+self.hz*60])
            ph_end = np.mean(self.ph_storage[self.j2[i]-self.hz*60:self.j2[i]])
            
            self.flux3o2_stor[i] = self.flux3o2[i] - (o2_start - o2_end) * self.B * self.mheight
            self.flux3ph_stor[i] = self.flux3ph[i] - (ph_start - ph_end) * self.B * self.mheight
            
            # Turbulent kinetic energy
            self.TKE3[i] = np.sqrt(np.mean(vx_prime**2) + np.mean(vy_prime**2) + 
                                   np.mean(vz_prime**2))
            
            # Wave velocity
            self.wv3[i] = np.sqrt(np.mean((self.vx[idx] - self.vxmean[i])**2) + 
                                 np.mean((self.vy[idx] - self.vymean[i])**2))
            
            # Cross-spectral density calculations
            # vx-vz cospectra
            freqxy, Cxy = signal.csd(vx_prime, vz_prime, fs=self.hz, 
                                     nperseg=self.window, scaling='density')
            df = freqxy[1] - freqxy[0]
            self.Flux_cpsd_xy[i] = np.sum(np.real(Cxy) * df)
            
            # Low frequency only
            freq_filt = freqxy < 1 / self.Td
            self.Flux_cpsd_xy_low[i] = np.sum(np.real(Cxy[freq_filt]) * df)
            
            # Friction velocity and roughness
            self.Ustar[i] = np.sqrt(np.abs(self.Flux_cpsd_xy_low[i])**0.5)
            if self.Ustar[i] > 0:
                self.Zzero[i] = self.mheight / np.exp((self.vmean[i] * 0.41) / 
                                                       (self.Ustar[i] + 1e-10))
                self.EddyDiff[i] = 0.41 * self.Ustar[i] * self.mheight
            
            # Drag coefficient (only for high flow)
            if self.vmean[i] >= self.flagwv:
                self.Cd[i] = np.abs(self.Flux_cpsd_xy_low[i]) / (self.vmean[i]**2)
                self.phi2[i] = self.phid[i]
            
            # O2 cospectra
            freqo2, Co2 = signal.csd(o2_prime, vz_prime, fs=self.hz, 
                                     nperseg=self.window, scaling='density')
            df = freqo2[1] - freqo2[0]
            self.Flux3cpsd[i] = np.sum(np.real(Co2) * df) * 3600
            
            freq_filt = freqo2 < 1 / self.Td
            self.Flux3cpsd_low[i] = np.sum(np.real(Co2[freq_filt]) * df) * 3600
            self.Flux3cpsd_low_stor[i] = self.Flux3cpsd_low[i] - \
                (o2_start - o2_end) * self.B * self.mheight
            
            # pH cospectra
            freqph, Cph = signal.csd(ph_prime, vz_prime, fs=self.hz, 
                                     nperseg=self.window, scaling='density')
            df = freqph[1] - freqph[0]
            self.Flux3cpsd_ph[i] = np.sum(np.real(Cph) * df) * 3600
            
            freq_filt = freqph < 1 / self.Td
            self.Flux3cpsd_ph_low[i] = np.sum(np.real(Cph[freq_filt]) * df) * 3600
            self.Flux3cpsd_ph_low_stor[i] = self.Flux3cpsd_ph_low[i] - \
                (ph_start - ph_end) * self.B * self.mheight
    
    def write_output(self):
        """Write results to output files."""
        print("Writing data to files...")
        
        # Prepare output data
        output_data = np.column_stack([
            self.timemean, self.vmean, self.o2mean, self.phmean, self.presmean, 
            self.Hs, self.TKE2, self.TKE3, self.flux2o2, self.flux3o2, 
            self.flux3o2_stor, self.Flux3cpsd_low, self.Flux3cpsd_low_stor,
            self.flux2ph, self.flux3ph, self.flux3ph_stor, 
            self.Flux3cpsd_ph_low, self.Flux3cpsd_ph_low_stor,
            self.SNRmean, self.correlationmean, self.vxmean, self.vymean, 
            self.vzmean, self.vxrot, self.vyrot, self.vzrot, self.thetad, 
            self.phid, self.phi2, self.horz, self.wv2, self.wv3, 
            self.Ustar, self.Cd, self.Zzero, self.EddyDiff, self.Flux_cpsd_xy_low
        ])
        
        # Select output file based on rotation
        if self.flagrotate == 0:
            outfile = f"outfile2_{Path(self.infile).name}"
        else:
            outfile = f"outfile2Rot_{Path(self.infile).name}"
        
        # Write to file
        header = ("timemean vmean o2mean phmean presmean Hs TKE2 TKE3 flux2o2 flux3o2 "
                 "flux3o2stor Flux3cpsdLow Flux3cpsdLowStor flux2ph flux3ph flux3phstor "
                 "Flux3cpsdphLow Flux3cpsdphLowStor SNRmean correlationmean vxmean vymean "
                 "vzmean vxrot vyrot vzrot theta phid phi2 horz wv2 wv3 Ustar Cd Zzero "
                 "EddyDiff Flux_cpsd_xy_low")
        
        # Use consistent format for all columns
        fmt_str = '%.6f'
        np.savetxt(outfile, output_data, fmt=fmt_str, header=header)
        
        print(f"Output written to: {outfile}")
    
    def process(self):
        """Run complete flux calculation."""
        try:
            self.load_data()
            self.find_burst_indices()
            self.convert_ph()
            self.calculate_burst_means()
            self.calculate_rotations()
            self.calculate_horizontal_angle()
            self.apply_timeshift()
            self.calculate_flux2()
            self.calculate_flux3()
            self.write_output()
            print("Done!")
            
        except Exception as e:
            print(f"Error: {e}", file=sys.stderr)
            sys.exit(1)


def main():
    """Main entry point for the script."""
    if len(sys.argv) < 3:
        print("Usage: python eddyflux_batch.py <data_file> <timestamp_file> [rotation_index]")
        print("  rotation_index: Required if using planar rotation (flagrotate=2)")
        sys.exit(1)
    
    infile = sys.argv[1]
    timestampfile = sys.argv[2]
    
    # Default to standard rotation
    flagrotate = 1
    rotation_index = None
    
    if len(sys.argv) >= 4:
        flagrotate = 2
        rotation_index = int(sys.argv[3])
    
    processor = EddyFluxProcessor(infile, timestampfile, flagrotate, rotation_index)
    processor.process()


if __name__ == "__main__":
    main()
