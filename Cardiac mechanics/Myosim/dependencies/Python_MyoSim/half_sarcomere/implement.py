# Functions for implementating half-sarcomere class
import numpy as np
import pandas as pd


def update_simulation(self, time_step, delta_hsl, hsl, y0, pf, cbf, calcium, n_array_length, cell_time,hs_params_new_list,set_data = 0):
    self.update_hs_props(hs_params_new_list)
    # time_step = time_step/1000.0
    self.Ca_conc = calcium
    self.temp_overlaps = 0.0
    self.hs_length = hsl
    self.myof.cb_force = cbf
    self.myof.pas_force = pf
    self.hs_force = self.myof.cb_force+self.myof.pas_force
    self.myof.y = y0
    if (np.abs(delta_hsl) > 0.0):
        self.myof.move_cb_distributions(delta_hsl)

    y_interp = self.myof.y
    self.myof.evolve_kinetics(time_step, self.Ca_conc, cell_time)
    # Assign int point's population vector to larger y vector
    y_pops = self.myof.y
    self.temp_overlaps = self.myof.n_overlap
    return self.temp_overlaps, y_interp, y_pops

def return_rates_fenics(self):
    return self.myof.return_fluxes(self.myof.y, self.Ca_conc)

def update_data_holder(self, dt, activation):

    # Update data struct for half-sarcomere
    self.data_buffer_index = self.data_buffer_index + 1
    self.hs_time = self.hs_time + dt
    self.hs_data.at[self.data_buffer_index, 'hs_time'] = self.hs_time
    self.hs_data.at[self.data_buffer_index, 'activation'] = activation
    self.hs_data.at[self.data_buffer_index, 'Ca_conc'] = self.Ca_conc
    self.hs_data.at[self.data_buffer_index, 'hs_length'] = self.hs_length
    self.hs_data.at[self.data_buffer_index, 'hs_force'] = self.hs_force
    self.hs_data.at[self.data_buffer_index, 'cb_force'] = self.myof.cb_force
    self.hs_data.at[self.data_buffer_index, 'pas_force'] = self.myof.pas_force

    if (self.myof.kinetic_scheme == '3state_with_SRX'):
        self.hs_data.at[self.data_buffer_index, 'M_OFF'] = \
            self.myof.y[0]
        self.hs_data.at[self.data_buffer_index, 'M_ON'] = \
            self.myof.y[1]
        self.hs_data.at[self.data_buffer_index, 'M_bound'] = \
            np.sum(self.myof.y[2 + np.arange(self.myof.no_of_x_bins)])
        self.hs_data.at[self.data_buffer_index, 'n_off'] = \
            np.sum(self.myof.y[-2])
        self.hs_data.at[self.data_buffer_index, 'n_on'] = \
            np.sum(self.myof.y[-1])

        # Update fluxes
        fluxes = self.myof.return_fluxes(self.myof.y, self.Ca_conc)
        self.hs_data.at[self.data_buffer_index, 'J1'] = fluxes['J1']
        self.hs_data.at[self.data_buffer_index, 'J2'] = fluxes['J2']
        self.hs_data.at[self.data_buffer_index, 'J3'] = np.sum(fluxes['J3'])
        self.hs_data.at[self.data_buffer_index, 'J4'] = np.sum(fluxes['J4'])
        self.hs_data.at[self.data_buffer_index, 'Jon'] = fluxes['Jon']
        self.hs_data.at[self.data_buffer_index, 'Joff'] = fluxes['Joff']

    if (self.myof.kinetic_scheme == '4state_with_SRX'):
        self.hs_data.at[self.data_buffer_index, 'M_OFF'] = \
            self.myof.y[0]
        self.hs_data.at[self.data_buffer_index, 'M_ON'] = \
            self.myof.y[1]
        self.hs_data.at[self.data_buffer_index, 'M_bound'] = \
            np.sum(self.myof.y[2 + np.arange(self.myof.no_of_x_bins)])
        self.hs_data.at[self.data_buffer_index, 'n_off'] = \
            np.sum(self.myof.y[-2])
        self.hs_data.at[self.data_buffer_index, 'n_on'] = \
            np.sum(self.myof.y[-1])

        # Update fluxes
        fluxes = self.myof.return_fluxes(self.myof.y, self.Ca_conc)
        self.hs_data.at[self.data_buffer_index, 'J1'] = fluxes['J1']
        self.hs_data.at[self.data_buffer_index, 'J2'] = fluxes['J2']
        self.hs_data.at[self.data_buffer_index, 'J3'] = np.sum(fluxes['J3'])
        self.hs_data.at[self.data_buffer_index, 'J4'] = np.sum(fluxes['J4'])
        self.hs_data.at[self.data_buffer_index, 'J5'] = np.sum(fluxes['J5'])
        self.hs_data.at[self.data_buffer_index, 'J6'] = np.sum(fluxes['J6'])
        self.hs_data.at[self.data_buffer_index, 'J7'] = np.sum(fluxes['J7'])
        self.hs_data.at[self.data_buffer_index, 'J8'] = np.sum(fluxes['J8'])
        self.hs_data.at[self.data_buffer_index, 'Jon'] = fluxes['Jon']
        self.hs_data.at[self.data_buffer_index, 'Joff'] = fluxes['Joff']        

    if (self.membr.kinetic_scheme == "Ten_Tusscher_2004"):
        # Ten Tusscher membrane voltage is in mV
        self.hs_data.at[self.data_buffer_index, 'membrane_voltage'] = \
            0.001*self.membr.y[0]

    self.hs_data.at[self.data_buffer_index, 'cb_number_density'] = self.cb_number_density
