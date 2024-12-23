## AdReDi - function packages
#
#
######################
## Python functions ##
######################
import numpy as np
import pandas as pd
import scipy as sc
from copy import deepcopy
from collections import defaultdict
import os
import math



##############
## CLASS -  ##
##############
class ShieldingGamma:

  def __init__(self,
               reactor_materials,
               source_params,
               geometry,
               mode,
               source_folder,
               nist_folder='NIST_att_coeff/'):
    """
    Class initialization

    Parameters
    ----------
    nist_path :
    """
    self.nist_names_map = {
        'Water':        'Water',
        'Air':          'Air',
        'Oxygen':       'O',
        'Nitrogen':     'N',
        'Potassium':    'K',
        'Beryllium':    'Be',
        'Boron':        'B',
        'Carbon':       'C',
        'Fluorine':     'F',
        'Iron':         'Fe',
        'Sodium':       'Na',
        'Uranium':      'U',
        'Zirconium':    'Zr',
        'Chromium':     'Cr',
        'Nickel':       'Ni',
        'Lead':         'Pb',
        'Bismuth':       'Bi',
        'Tungsten':     'W',
	'Concrete':	'Concrete',
	'Polyethylene':	'Polyethylene'
    }
    self.geometry = geometry
    self.mode = mode
    
    self.nist_folder = nist_folder
    self.reactor_materials = reactor_materials
    self.core_rho = 0
    self.materials_db = self.__materials_load_process__()

    if (source_folder is None) and (mode == 'full_analytic'):
      raise AttributeError('source_folder is required in full analytic mode')
    else:
      self.source_folder = source_folder
    self.source_db = source_params
    self.fluence_to_Gy_hr = 5.76e-7

    if self.mode == 'full_analytic':
      _ = self.source_initialization_full_analytic()
    elif self.mode == 'serpent_based':
      _ = self.source_initialization_serpent_based()
    elif self.mode == 'third_mode':
      _ = self.source_initialization_third_mode()
    else:
      raise AttributeError("mode must be either 'full_analytic', 'serpent_based' or 'third_mode'")

    _ = self.mu_rho_initialization()




  def __load_db__(self):
    nist_dict = defaultdict(pd.DataFrame)
    for nist_file in os.listdir(self.nist_folder):
      tmp_nist = pd.read_csv(
        f'{self.nist_folder}/{nist_file}',
        index_col=0
      )
      nist_dict[self.nist_names_map[nist_file[:-4]]] = tmp_nist
    return nist_dict





  def __materials_load_process__(self):

    # Load NIST database into a dictionary
    nist_db = self.__load_db__()

    # Process custom materials and elements
    if self.mode == 'full_analytic':
      salt_core = self.reactor_materials['salt_core']['rho'] * self.reactor_materials['salt_core']['vol'] / \
                  (self.reactor_materials['salt_core']['rho'] * self.reactor_materials['salt_core']['vol'] + \
                   self.reactor_materials['graphite_core']['rho'] * self.reactor_materials['graphite_core']['vol'])

      graphite_core = self.reactor_materials['graphite_core']['rho'] * self.reactor_materials['graphite_core']['vol'] / \
                      (self.reactor_materials['salt_core']['rho'] * self.reactor_materials['salt_core']['vol'] + \
                       self.reactor_materials['graphite_core']['rho'] * self.reactor_materials['graphite_core']['vol'])

      self.core_rho = self.reactor_materials['salt_core']['rho'] * salt_core + self.reactor_materials['graphite_core'][
        'rho'] * graphite_core
    elif self.mode == 'third_mode':
      i = 0
      vol_tot = 0
      for mat in self.reactor_materials.keys():
        if 'core' in self.reactor_materials[mat].keys():
          vol_tot += self.reactor_materials[mat]['vol']
      for mat in self.reactor_materials.keys():
        if 'core' in self.reactor_materials[mat].keys():
          self.core_rho += self.reactor_materials[mat]['rho'] * self.reactor_materials[mat]['vol']/vol_tot

    for elem in self.reactor_materials.keys():
      for ee, cc in enumerate(self.reactor_materials[elem]['components'].keys()):
        if ee == 0:
          fillv = 0
        else:
          fillv = None
        nist_db[elem] = nist_db[elem].add(self.reactor_materials[elem]['components'][cc] * nist_db[cc],
                                          fill_value=fillv).dropna()

    if self.mode == 'full_analytic':
      nist_db['core'] = salt_core * nist_db['salt_core'] + graphite_core * nist_db['C']
    elif self.mode == 'third_mode':
      nist_db['core'] = pd.read_csv(
        self.nist_folder + 'Air.csv',
        index_col=0,
      )

      nist_db['core'].loc[:, :] = 0
      e_scale = nist_db['core'].index.values
      for mat in self.reactor_materials.keys():
        if 'core' in self.reactor_materials[mat].keys():
          # in order to sum nist_db of different materials they are all interpolated on a single unified energy axis
          # taken from the air data
          nist_interp_mu = sc.interpolate.interp1d(nist_db[mat].index.values, nist_db[mat]['mu_rho (cm2/g)'])
          nist_interp_mut = sc.interpolate.interp1d(nist_db[mat].index.values, nist_db[mat]['mut_rho (cm2/g)'])

          nist_interp = np.zeros([len(e_scale), 2])
          nist_interp[:, 0] = nist_interp_mu(e_scale).reshape(-1)
          nist_interp[:, 1] = nist_interp_mut(e_scale).reshape(-1)
          # core homogenization by weighing each nist_db on the volume fraction
          nist_db['core'] += nist_interp * self.reactor_materials[mat]['vol']/vol_tot

    return nist_db

  def __merge_nested_dicts__(self, d1, d2):
    for key, value in d2.items():
      if key in d1 and isinstance(d1[key], dict) and isinstance(value, dict):
        self.__merge_nested_dicts__(d1[key], value)
      else:
        d1[key] = value




  def source_initialization_full_analytic(self):
    fname = self.source_db['fissile_material'] + '_gamma_spectrum.csv'
    df_load_s = pd.read_csv(
      f'{self.source_folder}/{fname}',
      index_col=0,
      sep=';'
    )

    if self.source_db['energy_res'] > 0:
      df_source = pd.DataFrame()
      df_source['Energy'] = np.linspace(df_load_s.index[0], df_load_s.index[-1], int(self.source_db['energy_res']))
      df_source['Energy Pdf'] = np.interp(df_source['Energy'], df_load_s.index, df_load_s['Gamma Fraction'])
      deltaE = df_source['Energy'].diff().iloc[-1]
      df_source['Gamma Spectrum'] = df_source['Energy Pdf'] * deltaE
      df_source['Energy upper'] = df_source['Energy'] / 1e6
      df_source['Energy lower'] = df_source['Energy'].shift(1) / 1e6
      df_source['Energy bin'] = df_source[['Energy lower', 'Energy upper']].mean(axis=1)
      df_source.drop(columns='Energy', index=0, inplace=True)
      df_source.set_index('Energy bin', drop=True, inplace=True)
    else:
      df_source = pd.DataFrame(index=df_load_s.index)
      deltaE = df_load_s.index.diff()
      df_source['Gamma Spectrum'] = df_load_s['Gamma Fraction'] * deltaE

    rv1 = self.source_db['spatial_res']
    rv2 = int(self.source_db['R_core'] + 1)
    drv = int(self.source_db['R_core'] / self.source_db['spatial_res'] + 1)
    semi_res = self.source_db['spatial_res'] / 2
    dr = np.linspace(rv1, rv2, drv) - semi_res

    df_bessel = pd.DataFrame(index=dr)
    norm_pow_dens = sc.special.jv(0, 2.4048 * dr / (self.source_db['R_core'] + self.source_db['extrap_r']))
    vol = 4 / 3 * np.pi * ((dr + semi_res) ** 3 - (dr - semi_res) ** 3)
    df_bessel['Norm Power'] = vol * norm_pow_dens

    gamma_source = np.outer(df_bessel['Norm Power'], df_source['Gamma Spectrum'])
    gamma_source = gamma_source / np.sum(gamma_source) * self.source_db['fission_rate'] * self.source_db[
      'photon_multiplicity']

    # if plot_flag:           #FIXME togliere i plot da qui, fare funzione specifica
    #   fig1 = mmp.mmp_plot(data=df_bessel, kind='line',
    #                       xlabel='Core Radius (cm)', ylabel='Bessel Norm Power', size='small')
    #   fig2 = mmp.mmp_plot(x=df_bessel.index, y=df_bessel["Norm Power"]/df_bessel['Norm Power'].sum(),
    #                       kind='bar', size='small',
    #                       xlabel='Core Radius (cm)', ylabel='Normalized Photon Production (1/cm3)')
    #   fig3 = mmp.mmp_plot(data=df_source, x=df_source.index, y='Gamma Spectrum',
    #                       kind='line', mode='lines+markers', log_x=True, log_y=True,
    #                       xlabel='Energy (MeV)', ylabel='Photon Energy Distribution (1/eV)',
    #                       size='small')

    # Save in self
    self.gamma_source = gamma_source
    self.energy_binning = df_source.index.values
    self.spatial_binning = df_bessel.index.values

    return


  def source_initialization_serpent_based(self):
    fname = self.source_db['energy_file_path']
    df_load_s = pd.read_csv(
        fname,
        index_col=0,
        sep=','
    )

    df_source = pd.DataFrame(index=df_load_s.index)
    df_source['Gamma Spectrum'] = df_load_s['Gamma Source']

    # Save in self
    self.gamma_source = df_source.values.reshape(1,-1)
    self.energy_binning = df_source.index.values
    self.spatial_binning = [self.source_db['source_radius']]

    return

  def source_initialization_third_mode(self):
    fname = self.source_folder + self.source_db['source_name']
    df_load_s = pd.read_csv(
      fname,
      index_col=0,
      sep=','
    )

    df_source = pd.DataFrame(index=df_load_s.index, )
    #de = np.append(0, df_load_s.index.values[1:] - df_load_s.index.values[0:-1])
    df_source['Gamma Spectrum'] = df_load_s['Gamma Source']/sum(df_load_s['Gamma Source'])*self.source_db['gamma_rate']

    # Save in self
    df_source.values.reshape(1, -1)
    self.spatial_binning = np.arange(0, self.source_db['spatial_binning_n'], dtype=float)
    dv = self.source_db['core_vol'] / self.source_db['spatial_binning_n']
    if self.geometry == 'sphere':
      r = (3 / 4 / np.pi * dv) ** (1 / 3)
      for i in range(0, self.source_db['spatial_binning_n']):
        self.spatial_binning[i] = r
        r = (3 / 4 / np.pi * dv + r ** 3) ** (1 / 3)
      self.energy_binning = df_source.index.values
      self.gamma_source = np.zeros([len(self.spatial_binning), len(self.energy_binning)])

      for i in range(0, self.source_db['spatial_binning_n']):
          self.gamma_source[i, :] = df_source.values.reshape(1, -1) / self.source_db['spatial_binning_n']
    else:
      r = (dv/self.source_db['H_core']/np.pi) ** (1/2)
      for i in range(0, self.source_db['spatial_binning_n']):
        self.spatial_binning[i] = r
        r = (dv/self.source_db['H_core']/np.pi + r ** 2) ** (1/2)
      self.energy_binning = df_source.index.values
      self.gamma_source = np.zeros([len(self.spatial_binning), len(self.energy_binning)])

      for i in range(0, self.source_db['spatial_binning_n']):
        self.gamma_source[i, :] = df_source.values.reshape(1, -1) / self.source_db['spatial_binning_n']
    #fig = go.Figure()
    #fig.add_trace(
    #  go.Scatter(x=df_load_s.index.values, y=self.gamma_source[0, :],
    #             line=dict(color='black', dash='solid', width=3.5), mode='lines')
    #)
    #fig.update_layout(title='', xaxis=dict(title='energy (MeV)'), yaxis=dict(title='Activity (Bq)'))
    #fig.show()

    return

  def mu_rho_initialization(self):

    # Log-Log interpolationg function
    def mu_rho(energy, df_ref, method='log'):
      idx_name = df_ref.index.name
      if method == 'log':
        df_ref_new = np.log(df_ref.reset_index()).set_index(idx_name)
        energy_new = np.log(energy)
      elif method == 'linear':
        df_ref_new = df_ref.copy()
        energy_new = energy.copy()
      values = pd.concat([df_ref_new, pd.DataFrame(index=np.array(energy_new))]).sort_index().interpolate(method='index').loc[energy_new].reset_index()
      if method == 'log':
        res = np.exp(values).set_index('index')
      elif method == 'linear':
        res = values.set_index('index')
      return res

    mu_rho_db = defaultdict(np.array)
    mut_rho_db = defaultdict(np.array)
    for elem in list(self.reactor_materials.keys()):
      mu_rho_db[elem] = mu_rho(self.energy_binning, self.materials_db[elem])['mu_rho (cm2/g)'].values
    mut_rho_db['Air'] = mu_rho(self.energy_binning, self.materials_db['Air'])['mut_rho (cm2/g)'].values
    if self.mode == 'full_analytic' or self.mode == 'third_mode':
      mu_rho_db['core'] = mu_rho(self.energy_binning, self.materials_db['core'])['mu_rho (cm2/g)'].values

    # Save in self
    self.mu_rho_db = mu_rho_db
    self.mut_rho_db = mut_rho_db

    return




  def dose_calc(self,
                case_geom,
                target_distance):

    mat_dose = np.zeros_like(self.gamma_source)

    # Update materials dictionary
    case_materials_db = deepcopy(self.reactor_materials)
    self.__merge_nested_dicts__(case_materials_db, case_geom)  #######

    # Calculate mat_dose
    for eni, en in enumerate(self.energy_binning):

      materials_att = 0
      for mr in self.reactor_materials.keys():
        if 'core' not in self.reactor_materials[mr].keys() and self.mode == 'third_mode':
            materials_att += self.mu_rho_db[mr][eni] * case_materials_db[mr]['rho'] * case_materials_db[mr]['thickness']
        elif mr not in ['salt_core', 'graphite_core', 'core', 'Air'] and self.mode != 'third_mode':
            materials_att += self.mu_rho_db[mr][eni] * case_materials_db[mr]['rho'] * case_materials_db[mr][
                'thickness']
      if math.isnan(materials_att):
        materials_att = 0

      for ri, r in enumerate(self.spatial_binning):
        if self.mode == 'serpent_based':
          core_att = 0
        else:
          core_att = self.mu_rho_db['core'][eni] * self.core_rho * (self.source_db['R_core'] - r)
        exponent = - (core_att + materials_att)
        if self.geometry == 'sphere':
          mat_dose[ri, eni] = self.fluence_to_Gy_hr * self.gamma_source[ri, eni] * en * np.exp(exponent) / \
                              (4 * np.pi * target_distance ** 2) * self.mut_rho_db['Air'][eni]
        else:
          mat_dose[ri, eni] = self.fluence_to_Gy_hr * self.gamma_source[ri, eni] * en * np.exp(exponent)/ \
          (2 * np.pi * (target_distance * self.source_db['H_core'])) * self.mut_rho_db['Air'][eni]


    df_dose = pd.DataFrame(index=self.spatial_binning, data=mat_dose, columns=self.energy_binning)
    self.df_dose = df_dose

    tot_dose = df_dose.sum().sum()#*1e6
    self.tot_dose = tot_dose


    return self.tot_dose  #in Gy/h
