"""Here are the action items that can be taken by the user to control
the spectrometer and the stage.

"""

import numpy as np
import time
from tkinter import messagebox

from outsource import outsource
from ui_plot import update_plot, ui_plot_graph, get_smooth_data, show_legend
from prog_bar import prog_window
from suzy_constants import ORIGIN_POS, LOAD_POS, WAVENUMBER_RANGE, SETT_OPTIMIZE
from suzy_constants import CADDY_SIZE, OPT_THRESHOLD, OPT_SPAN, OPT_DELTA
from suzy_constants import OPT_SPOT, OPT_ERR, DEBUG_SERIAL, STAGE_CMD_DELAY
from suzy_constants import DT_SHOW_RAW, CELL_RMAP
from master import get_dataSetID, add_dataFile, add_dataSetEntry
from fileIO import make_dataFile_obj
from dialogs import entry_prompt


X, Y, Z, REL, ABS = 'X', 'Y', 'Z', 'REL', 'ABS'

# check that hardware connection is possible
try:
    from logan import scanSierra
    NO_HW = False

except Exception as e:
    print('Error: ', e, '\nHardware not available.')
    NO_HW = True

    
class controls(object):
    """This class contains all the action items for the UI."""

    def __init__(self, ui_parent, origin_pos=ORIGIN_POS,
                 load_pos=LOAD_POS, xr=WAVENUMBER_RANGE,
                 p_sett=SETT_OPTIMIZE, caddy_size=CADDY_SIZE,
                 opt_thresh=OPT_THRESHOLD, opt_span=OPT_SPAN,
                 opt_delta=OPT_DELTA, opt_spot=OPT_SPOT,
                 stage_cmd_delay=STAGE_CMD_DELAY, opt_err=OPT_ERR):
        """Setup the controls structure to work with the spectrometer and the
        stage.

        Args:
            ui_parent: the pointer to the controls UI class

        """
        # set up calibrated coordinates for A1 caddy element
        self.origin = np.array(origin_pos)
        self.load = np.array(load_pos)
        # set up wavenumber array
        self.xr = np.array(xr)

        self.p_sett = p_sett
        self.caddy_size = caddy_size
        self.opt_thresh = opt_thresh
        self.opt_span = opt_span
        self.opt_delta = opt_delta
        self.opt_spot = opt_spot
        self.opt_err = opt_err
        self.stage_cmd_delay = stage_cmd_delay
        self.z_reset = True

        # tentative resource names
        self.ports = ['COM4', 'COM3']

        self.ui = ui_parent
        self.curve = None

        # labels
        self.set_lbl = None
        self.spec_lbl = None
        self.set_tlbl = None
        self.spec_id = 0

        
    def connect(self, ui, on=True):
        """This function connects to the hardware.

        Args:
            ui: the pointer to the main window UI
    
            on: the flag to either connect (True) or disconnect (False).
                Defaults to True

        Returns: 
            A list [boolean,boolean] corresponding to the
            success in connection to the spectrometer and stage
            respectively.
        """

        if NO_HW:
            return False, False
        
        pbar = prog_window(ui, ttl='CONNECTING TO HARDWARE' if on \
                           else 'DISCONNECTING FROM THE HARDWARE',
                           allow_abort=False)
        
        if on:
            # connect to stage and spectrometer
            self.logan = scanSierra(self, self.origin, self.ports,
                                    DEBUG_SERIAL, prog=pbar)
            
            # update resource names
            self.ports = self.logan.ports

            # devices status: [stage, spectrom.]
            self.enabled = self.logan.enabled

            # update the current position label
            self.update_pos_lbl()

        else:
            pbar.start(3)
            pbar.set_status('Resetting the spectrometer ...')
            self.reset_device()
            pbar.update(1)
            
            if hasattr(self, 'logan'):
                self.logan.close(pbar=pbar)
                
            self.enabled = [False, False]
            
        pbar.destroy()
        return self.enabled[0], self.enabled[1]


    def setup_device(self):
        """Setup the device according to the user settings. 

        TURN LASER ON
        """
        if not self.enabled[0]:
            return
        
        # clear the plot figure
        self.curve = None

        # update the user settings
        self.sett = self.ui.get_settings()

        # setup spectrometer
        self.logan.settings(BAS='OFF',
                            RAS=self.sett[0],
                            REF=self.sett[1],
                            AVE=int(self.sett[2]),
                            PLV=int(self.sett[3]),
                            INT=float(self.sett[4]),
                            PWR='ON', SHU='OPEN')
        

    def setup_optimize(self):
        """Setup the device according to 'producer' settings
        for optimization routines
        """
        if not self.enabled[0]:
            return

        # set up the device for optimization measurements
        self.logan.settings(**self.p_sett)

        
    def reset_device(self):
        """This function resets the device and turns off the laser.

        TURN LASER OFF
        """
        if not hasattr(self, 'enabled') or not self.enabled:
            return

        self.logan.settings(SHU='CLOSE', PWR='OFF', RAS='OFF')

        
    def read_spectrum(self, ui, run_setup=False):
        """This function triggers a NEW ACQUISITION from the spectrometer.

        Args:
            ui: the pointer to the MAIN WINDOW ui.

            run_setup: set True to reset the device, False to
                ignore. Defaults to False.

        Returns:
            A tuple with the (x, y) for the new acquired spectrum.

        """

        if not self.enabled[0]:
            yr = np.zeros(1801) 
            return self.xr, yr

        if run_setup:
            pbar = prog_window(ui, ttl='ACQUIRING DATA', allow_abort=False)
            pbar.start(3)
            pbar.set_status('Setting up hardware...')
            self.setup_device()
        else:
            pbar = None
        
        if pbar is not None:
            pbar.update(1)
            pbar.set_status('Acquiring data...')

        yr = self.logan.acquire()

        if run_setup:
            if pbar is not None:
                pbar.update(2)
                pbar.set_status('Shutting down the hardware ...')

            self.reset_device()

            if pbar is not None:
                pbar.destroy()
                
        return self.xr, yr


    def tree_insert_spectrum(self, ui, yr):
        q = entry_prompt(ui,
                         ['Measurement set:', 'Spectrum label:'],
                         [self.set_lbl if self.set_lbl is not None else '',
                          self.spec_lbl if self.spec_lbl is not None else ''],
                         'Enter the acquisition information.',
                         field_len=30)

        if q.result is not None:
            mlbl, slbl = q.result

            if self.set_lbl is None or mlbl != self.set_lbl:
                # start a new measurement set
                self.set_lbl = mlbl if len(mlbl) else 'SNRI'
                self.spec_id = 0
                
                # get the next available entry id
                set_id, self.set_tlbl = get_dataSetID(ui)
                
                # make a top-level entry for this data set:
                add_dataSetEntry(ui, set_id, self.set_tlbl, self.set_lbl)

            else:
                # add to an existing set
                self.spec_id += 1

            self.spec_lbl = slbl if len(slbl) else 'spectrum'
            # prompt the user to save the spectrum
            add_dataFile(ui, self.set_tlbl, self.spec_id,
                         self.set_tlbl+'::spectrum'+str(self.spec_id),
                         make_dataFile_obj(self.xr, yr, self.spec_lbl,
                                           master=ui))


    def reset_origin(self):
        try:
            x = float(self.logan.stage.get_position('X'))
            y = float(self.logan.stage.get_position('Y'))
            z = float(self.logan.stage.get_position('Z'))

            # update the cell psotions
            self.logan.set_cell_pos([x, y-2*12.75, z])
            
        except:
            pass

        
    def update_pos_lbl(self):
        """Show the position label in the UI."""

        try:
            x = float(self.logan.stage.get_position('X'))
            y = float(self.logan.stage.get_position('Y'))
            z = float(self.logan.stage.get_position('Z'))
        
            self.ui.pos_lbl.set('x: {:.2f}    y: {:.2f}    z: {:.2f}'\
                                .format(x, y, z))
        except:
            pass

        
    def move(self, direction):
        """This function will make the move.

        Args:
            direction: a string indicating the direction to move in. Can be
                '+X', '-Z', '+Y', ...

        """
        if not self.enabled[1]:
            return
        
        # update the user settings
        self.sett = self.ui.get_settings()

        sign = 1 if direction[0] == '+' else -1
        axis = direction[1]

        # calculate the displacement requested
        dx = sign*float(self.sett[7 if axis == 'Z' else 6])

        if axis != 'Z' \
           or float(self.logan.stage.get_position('Z'))+dx < 10:
            self.logan.stage.move(axis, dx, REL, False)
            self.z_reset = False

        # update the current position label
        self.update_pos_lbl()
        
            
    def goto_load(self):
        """This function moves the stage to the LOAD positions."""
        if not self.enabled[1]:
            return
        
        self.logan.stage.move('X', self.load[0], ABS, False)

        time.sleep(self.stage_cmd_delay)
        self.logan.stage.move('Y', self.load[1], ABS, False)

        if self.z_reset:
            time.sleep(self.stage_cmd_delay)
            self.logan.stage.move('Z', self.load[2], ABS, False)
        
        # update the current position label
        self.update_pos_lbl()


    def goto_cell(self, loc='C1', z_dim = False):
        """This function moves the stage to the specified location at the
        center of the cell.

        Args:
            loc: the string indicating the cell ID to go to.
            
            z_dim: a flag to move or not the z dimention. Default False.
        """
        if not self.enabled[1]:
            return

        if self.z_reset:
            self.logan.goto_caddy(loc, False, True)
        else:
            self.logan.goto_caddy(loc, False, False)
            
        self.curr_cell = loc

        # update the current position label
        self.update_pos_lbl()

        
    def prep_hw_for_search(self, ui):
        # get the current settings
        self.sett = self.ui.get_settings()

        pbar = prog_window(ui, ttl='LOCATING CHIPS')
        pbar.set_status('Setting up the hardware')

        # prepare the hardware for search
        self.setup_device()
        self.setup_optimize()
        self.logan.stage.set_accel(slow=True)

        return pbar


    def reset_hw_from_search(self, pbar):
        # reset the device
        pbar.set_status('Shutting down the hardware')
        self.logan.stage.set_accel(slow=False)
        self.reset_device()
        
        # update the current position label
        self.update_pos_lbl()

        # keep current z-position
        self.z_reset = False

        # close the progress bar
        pbar.destroy()

    
    def run_find_chips(self, ui):
        """This function will locate the chips within a current cell.

        Args:
            ui: the pointer to the MAIN WINDOW UI.

        """
        if (not self.enabled[0]) or (not self.enabled[1]):
            return

        pbar = self.prep_hw_for_search(ui)                
        ydata = self.logan.find_samples(self.curr_cell, int(self.sett[5]),
                                        self.caddy_size, self.opt_spot,
                                        self.opt_err, self.opt_thresh,
                                        self.opt_span, self.opt_delta,
                                        pbar=pbar)
        self.reset_hw_from_search(pbar)

        return np.array(ydata)
    
    
    def run_auto_focus(self, ui):
        """This function will run an autofocus routine (only in Z).

        Args:
            ui: the pointer to the MAIN WINDOW UI.

        """
        if (not self.enabled[0]) or (not self.enabled[1]):
            return

        pbar = prog_window(ui, ttl='RUNNING AUTOFOCUS', allow_abort=False)

        # prepare the hardware for search
        pbar.start(5)
        pbar.set_status('Setting up the devices')
        
        self.setup_device()
        self.setup_optimize()
        self.logan.stage.set_accel(slow=True)

        pbar.update(1)
        z_best = self.logan.autofocus(ui, pbar, 1)

        pbar.update(4)
        pbar.set_status('Shutting down the hardware')
        self.reset_device()
        self.logan.stage.set_accel(slow=False)
        
        # keep current z-position
        self.z_reset = False
        
        # update the current position label
        self.update_pos_lbl()

        # close the progress bar
        pbar.destroy()
        
        return z_best

    
    def run_scan_caddy(self, ui):
        """This function will scan the whole caddy.

        Args:
            ui: the pointer to the MAIN WINDOW UI.

        """

        pbar = self.prep_hw_for_search(ui)
        pbar.start(len(self.sett[8]))

        p, lbls, ydata = [], [], []
        for i,cell in enumerate(self.sett[8]):
            cell_lbl = CELL_RMAP[cell]
            
            pbar.set_status('Scanning cell ' + cell_lbl)
            spectra = self.logan.find_samples(cell, int(self.sett[5]),
                                               self.caddy_size, self.opt_spot,
                                               self.opt_err, self.opt_thresh,
                                               self.opt_span, self.opt_delta,
                                               pbar=pbar)
            # accumulate the spectra
            ydata += spectra

            # accumulate the labels
            [lbls.append(cell_lbl) for each in range(len(spectra))]
            
            if pbar.is_aborted():
                break
            else:
                pbar.update(i+1)

        self.reset_hw_from_search(pbar)
        return np.array(ydata), lbls


    def save_caddy_spectra(self, res, ui):
        """This function will put the caddy data into the data tree.

        Args:
            res: a tuple with e NumPy array of dat and the list of labels

            ui: the pointer to the main UI

        """

        self.ui.disable_buttons('scan_btn', False)

        # split the result into data and labels
        data, lbls = res

        if not len(lbls):
            return

        # split up the labels
        lbl0, idx = lbls[0], [0]
        for i in range(1, len(lbls)):
            if lbls[i] != lbl0:
                idx.append(i)
                lbl0 = lbls[i]
        idx.append(i+1)

        # push the data to the tree
        for i in range(1, len(idx)):
            # create a new data node
            self.init_data_tree_node(ui, False, lbls[idx[i-1]])
            
            # populate the data tree
            self.push_data2tree(data[idx[i-1]:idx[i]], ui)

        messagebox.showinfo('SCAN DONE',
                            'Caddy scan is complete',
                            parent=ui)


            
    def init_data_tree_node(self, ui, prompt_user=True, mlbl=None):
        """This function will create a new higher-level node entry in the data
        tree.

        Args:
            ui: the pointer to the main UI

            prompt_user: offer a user a choice to name the node. Defaults 
                to True.

        Returns:
            True on sucess, False on user cancel.

        """
        
        if prompt_user:
            q = entry_prompt(ui,
                             ['Measurement set:'],
                             [self.set_lbl if self.set_lbl is not None else ''],
                             'Enter the acquisition information.',
                             field_len=30)

            # bail is user also bails
            if q.result is None:
                return False

            mlbl = q.result[0]

        if prompt_user == False \
           or self.set_lbl is None or mlbl != self.set_lbl:
            # start a new measurement set
            self.set_lbl = mlbl if len(mlbl) else 'SNRI'
            self.spec_id = 0
                
            # get the next available entry id
            set_id, self.set_tlbl = get_dataSetID(ui)
                
            # make a top-level entry for this data set:
            add_dataSetEntry(ui, set_id, self.set_tlbl, self.set_lbl)

        else:
            # add to an existing set
            self.spec_id += 1

        return True
    

    def push_data2tree(self, res, ui):
        """This function will populate the data tree.

        Args:
            res: a NumPy array with the data that makes up a data set.

            ui: a pointer to the main UI.

        """
        
        axis, count, raw = [], [], []

        for i in range(res.shape[0]):
            spec_lbl = 'chip '+str(i+1)

            # make a dataFile record out of the data
            rec = make_dataFile_obj(self.xr, res[i,:], spec_lbl,
                                    master=ui)
            add_dataFile(ui, self.set_tlbl, self.spec_id,
                         self.set_tlbl+'::spectrum'+str(self.spec_id),
                         rec)
            self.spec_id += 1

            # also, update the averages
            if len(axis) and len(axis) == len(rec.axis()):
                counts = np.vstack((counts, rec.counts()))
                raw = np.vstack((raw, rec.counts(data_type='raw')))
                    
            elif not len(axis):
                axis = rec.axis()
                counts = rec.counts()
                raw = rec.counts(data_type='raw')
                
        # also insert the average
        if res.shape[0] > 1:
            add_dataFile(ui, self.set_tlbl, self.spec_id,
                         self.set_tlbl+'::average'+str(self.spec_id),
                         make_dataFile_obj(axis,
                                           counts.T.mean(axis=1),
                                           'average',
                                           std=counts.T.std(axis=1),
                                           raw=raw.T.mean(axis=1),
                                           raw_std=raw.T.std(axis=1),
                                           master=ui))
            self.spec_id += 1
        

    def save_chips_spectra(self, res, ui):
        if self.init_data_tree_node(ui, prompt_user=True, mlbl=None):
            # populate the data tree
            self.push_data2tree(res, ui)

            # plot the results
            self.plot_spectrum(res, ui, tag='ui_find', multi=True)

        
    def find_chips(self):
        """This is a wrapper function to locate the chips in the current cell.

        """
        
        self.ui.disable_buttons('ui_find', True)
        outsource(self.run_find_chips, self.save_chips_spectra, self.ui.parent)

        
    def auto_focus(self):
        """This is a wrapper function to run the autofocus.

        """
        self.ui.disable_buttons('ui_af', True)
        outsource(self.run_auto_focus,
                  lambda res, ui: self.ui.disable_buttons('ui_af', False),
                  self.ui.parent)
        

    def scan_caddy(self):
        """This is a wrapper function to scan the whole caddy.

        """        
        self.ui.disable_buttons('scan_btn', True)

        outsource(self.run_scan_caddy, self.save_caddy_spectra,
                  self.ui.parent)


    def live_view(self):
        """Run the live view, constantly scanning the spectrum."""
        outsource(self.run_live, lambda res, ui: True, self.ui.parent)
    
        
    def run_live(self, ui):
        """This function actually reads the spectra.

        Args:
            ui: the pointer to the MAIN WINDOW UI.

        """

        
        self.setup_device()
        self.logan.settings(REF='OFF')
        
        while ui.controls.live:
            self.plot_spectrum(self.read_spectrum(ui), ui)

        self.reset_device()


    def acquire(self):
        """This function takes one acquisition from the device."""

        self.sett = self.ui.get_settings()        
        self.ui.disable_buttons('ui_acq', True)
        self.get_spectrum('ui_acq')

        
    def get_spectrum(self, tag=None):
        """This function pulls the spectrum from the device.

        Args
            tag: the UI element tag, optionally, to reset the button upon
                read completion. Defaults to None.

        """
        
        outsource(self.read_spectrum, self.plot_spectrum, self.ui.parent,
                  cmd_args=(True,), cleanup_args=(tag, True))

        
    def plot_spectrum(self, res, ui, tag=None, push_data=False,
                      multi=False):
        """This function plots the spectrum on the UI.

        Args:
            res: a tuple with the (x, y) values for the spectra

            ui: the pointer to the MAIN WINDOW UI

            tag: the optional tag for the ui element to reset. 
                Defaults to None.

            push_data: a flag to signal offering user to push the data
                to the data tree. Defaults to false

            multi: the plot type, either False for one curve, or
                True for multi-curve figure. Defaults to False.

        """
        # choose the data type to show
        dt = 'raw' if ui.data_selBox.current() == DT_SHOW_RAW else 'sig'
                
        # calculate the baseline data
        if not multi:
            rec = make_dataFile_obj(*res, 'tmp', master=ui)
            # check that the data was loaded correctly
            if rec is None:
                return

            xr, yr = rec.axis(), get_smooth_data(ui, rec, dt)
        
            if self.curve is None:
                ui.mpl_ax.cla()
                update_plot(ui)

                self.curve = ui_plot_graph(ui, xr, yr)[0]
            else:
                self.curve.set_ydata(yr)
                ymin, ymax = min(yr)-100, max(yr)+100
                ui.mpl_ax.set_ylim([ymin, ymax])

        else:
            xvals, Y = self.xr, res

            # clear the plot
            ui.mpl_ax.cla()
            update_plot(ui)
            
            # plot all the graphs
            for i,yvals in enumerate(Y):
                rec = make_dataFile_obj(xvals, yvals, 'tmp', master=ui)
                xr, yr = rec.axis(), get_smooth_data(ui, rec, dt)               
                ui_plot_graph(ui, xr, yr, label='chip '+str(i+1),
                              lc=ui.color_names[i%len(ui.color_names)])

            # put up the legend
            if len(Y):
                show_legend(ui)
                
                
        # update the figure
        ui.figCanvas.show()

        # optionally, reset the button that called the action
        if tag is not None:
            ui.controls.disable_buttons(tag, False)
        

        if push_data:
            try:
                self.tree_insert_spectrum(ui, yr)
            except Exception as e:
                print(e)
