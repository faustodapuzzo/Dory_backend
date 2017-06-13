"""This class administers the UI for the spectrometer/stage controls.

"""

import tkinter as tk
from tkinter import ttk
from tkinter import N, W, E, S
from tkinter import messagebox

from controls_ai import controls
from suzy_constants import CELL_NAMES, CELL_IDS, AVE_OPTIONS, FEATURE_SET4
from suzy_constants import CELL_MAP
from outsource import outsource

class controls_ui(object):
    """This class adds the UI widgets to allow user control of the
    spectrometer.

    """
    def __init__(self, parent, frame):
        """Place the UI elements for the spectrometer control.

        Args:
            parent: the reference to the main UI

            frame: the main frame where the UI controls should go.

        """

        # save the state variables
        self.live = False
        self.parent = parent
        self.mf = frame

        # prepare the structure to handle the action items
        self.ai = controls(self)

        # place the controls
        self.placeControls()

        # get the initial settings
        self.get_settings()

        
    def get_settings(self):
        """This functions grabs the settings for the acquisition from the UI.
    
        """
        cells = []
        for col in CELL_NAMES:
            for row in CELL_IDS:
                if self.cell_var[col+row].get():
                    cells.append(CELL_MAP[col+row])
        
        self.sett = [self.ui_rstr.get(), 'OFF',
                     self.ui_ave.get(), self.ui_pwr.get(),
                     self.ui_int.get(), self.ui_nchips.get(),
                     self.ui_hstep.get(), self.ui_vstep.get(),
                     cells]
        
        return self.sett
    
        
    def layout_CtlCombo(self, frame_tag, var_tag, values):
        """This function places a combo box for the device control.

        Args:
            frame_tag: a string for identifying the frame to place 
                the combo box.

            var_tag: a string to identify the new combo box

            values: the values for the combo box

        Returns:
            The pointer to the new combo box

        """

        self.ctl[var_tag] = ttk.Combobox(self.ctl[frame_tag],
                                         state='readonly',
                                         values=values,
                                         width=8)
        self.ctl[var_tag].current(0)
        self.ctl[var_tag].bind('<<ComboboxSelected>>',
                               lambda evt: evt.widget.master.focus_set())

        return self.ctl[var_tag]


    def make_controlFrames(self):
        """Make all the frames for the control UI."""
        self.ctl = {}

        self.ctl['sett_fr'] = tk.Frame(self.mf, relief=tk.GROOVE,
                                       borderwidth=1)
        self.ctl['spec_fr'] = tk.Frame(self.mf)
        self.ctl['navi_fr'] = tk.Frame(self.mf)
        self.ctl['auto_fr'] = tk.Frame(self.mf)
        self.ctl['cell_fr'] = tk.Frame(self.mf)

        ttk.Separator(self.mf, orient=tk.HORIZONTAL)\
           .grid(column=0, row=1, padx=5, pady=5, columnspan=3, sticky=(W,E))

        self.ctl['sett_fr'].grid(column=0, row=2, padx=7, pady=2,
                                 sticky='nsew')
        self.ctl['spec_fr'].grid(column=1, columnspan=2, row=2, padx=7, pady=2,
                                 sticky='nsew')

        ttk.Separator(self.mf, orient=tk.HORIZONTAL)\
           .grid(column=0, row=3, padx=5, pady=5, columnspan=3, sticky=(W,E))

        self.ctl['navi_fr'].grid(column=0, columnspan=3, row=4, padx=7, pady=2,
                                 sticky='nsew')

        ttk.Separator(self.mf, orient=tk.HORIZONTAL)\
           .grid(column=0, row=5, padx=5, pady=5, columnspan=3, sticky=(W,E))

        self.ctl['auto_fr'].grid(column=0, columnspan=3, row=6, padx=7, pady=2,
                                 sticky='nsew')

        ttk.Separator(self.mf, orient=tk.HORIZONTAL)\
           .grid(column=0, row=7, padx=5, pady=5, columnspan=3, sticky=(W,E))

        self.ctl['cell_fr'].grid(column=0, columnspan=3, row=8, padx=7, pady=2,
                                 sticky='nsew')



    def add_CtlSettFrame(self):
        """Draw the controls to adjust the spectrometer settings."""
        
        ttk.Label(self.ctl['sett_fr'], text='Raster')\
            .grid(column=0, row=0, padx=2, pady=5, sticky='e')
        ttk.Label(self.ctl['sett_fr'], text='Reference', state='disabled')\
            .grid(column=0, row=1, padx=2, pady=5, sticky='e')
        ttk.Label(self.ctl['sett_fr'], text='Average')\
            .grid(column=0, row=2, padx=2, pady=5, sticky='e')

        self.ui_rstr = self.layout_CtlCombo('sett_fr', 'RASTR', ['ON', 'OFF'])
        self.ui_rstr.grid(column=1, row=0, padx=2, pady=5, sticky='e')

        self.ui_ref = self.layout_CtlCombo('sett_fr', 'REF', ['ON', 'OFF'])
        self.ui_ref['state'] = 'disabled'
        self.ui_ref.current(1)
        self.ui_ref.grid(column=1, row=1, padx=2, pady=5, sticky='e')

        self.ui_ave = self.layout_CtlCombo('sett_fr', 'AVE', AVE_OPTIONS)
        self.ui_ave.grid(column=1, row=2, padx=2, pady=5, sticky='e')


    def val_num(self, str_value, vmin, vmax):
        try:
            num = float(str_value)
            return num >= vmin and num <= vmax
        except:
            pass

        return False


    def live_view(self):
        """This functions engages the live view button."""
        # switch the state of the button

        self.ui_lv['text'] = 'LIVE VIEW' if self.live else 'STOP'
        self.disable_buttons('ui_lv', engage=not self.live)
        self.live = not self.live

        if self.live:
            self.ai.live_view()


    def connect(self, on=True):
        self.disable_buttons('ui_conn', True)

        outsource(self.ai.connect,
                  lambda res, ui: self.attempt_connect(res, on),
                  self.parent,
                  cmd_args=(on,))

        
    def attempt_connect(self, res, on=True):
        """This function initializes the connection to the instruments."""

        # attempt the connection to the hardware
        spec, stage = res
        spec_btn_state = 'normal' if spec else 'disabled'
        stage_btn_state = 'normal' if stage else 'disabled'
        
        if stage or spec:
            # reset the button state
            self.disable_buttons('ui_conn', False)

            # rename the button
            self.ui_conn['text'] = 'DISCONNECT' if on \
                                   else 'CONNECT TO SPECTROMETER'
            self.ui_conn['relief'] = 'sunken' if on else 'raised'
            self.ui_conn['command'] = lambda: self.connect(not on)

        else:
            # reset the button view
            self.ui_conn['text'] = 'CONNECT TO SPECTROMETER'
            self.ui_conn['relief'] = 'raised'
            self.ui_conn['command'] = lambda: self.connect(True)
            self.ui_conn['background'] = self.ui_mv[0]['background']


        if on and not stage:
            messagebox.showerror('Connection ERROR',
                                 'No stage found. Check the COM ports',
                                 parent=self.parent)

        if on and not spec:
            messagebox.showerror('Connection ERROR',
                                 'No spectrometer found. Check the COM ports',
                                 parent=self.parent)

        # enabke/disable the buttons
        self.ui_lv['state'] = spec_btn_state
        self.ui_acq['state'] = spec_btn_state

        for q in self.cell_btn.keys():
            self.cell_btn[q]['state'] = stage_btn_state

        for each in self.ui_mv:
            each['state'] = stage_btn_state

        for each in [self.ui_af, self.ui_find, self.gt_load,
                     self.gt_a1, self.scan_btn]:
            each['state'] = stage_btn_state



    def disable_buttons(self, tag, engage=True):
        """This function disables colliding actions.

        Args:
            tag: the tag for the button that is being activated.

        """

        # avoid any surprises
        if not hasattr(self, tag):
            return

        # disable confilicting buttons
        for each in [self.ui_lv, self.ui_acq, self.scan_btn,
                     self.ui_af, self.ui_find, self.ui_conn]:
            if each != getattr(self, tag):
                each['state'] = 'disabled' if engage else 'normal'
        
        # change the background
        getattr(self, tag)['relief'] = 'sunken' if engage else 'raised'
        getattr(self, tag)['background'] = 'red' if engage\
                                           else self.ui_mv[0]['background']

                
    def add_CtlSpecFrame(self):
        """Draw the spectrometer controls."""
        
        self.ui_lv = tk.Button(self.ctl['spec_fr'], text='LIVE VIEW',
                               state='disabled', command=self.live_view)
        self.ui_lv.grid(column=0, row=0, sticky='nsew', padx=5, pady=5)
        self.ui_acq = tk.Button(self.ctl['spec_fr'], text='ACQUIRE',
                                state='disabled', command=self.ai.acquire)
        self.ui_acq.grid(column=1, row=0, sticky='nsew', padx=5, pady=5)

        ttk.Label(self.ctl['spec_fr'], text='Power level')\
            .grid(column=0, row=1, padx=2, pady=5, sticky='e')
        ttk.Label(self.ctl['spec_fr'], text='Integration (s)')\
            .grid(column=0, row=2, padx=2, pady=5, sticky='e')

        self.ui_pwr = self.layout_CtlCombo('spec_fr', 'PWR',
                                           list(map(str, range(1,16))))
        self.ui_pwr.grid(column=1, row=1, padx=2, pady=5, sticky='w')

        vcmd = self.ctl['spec_fr']\
                   .register(lambda s: self.val_num(s, 0.1, 1000))
        self.ui_int = ttk.Entry(self.ctl['spec_fr'],
                                validate='focus',
                                validatecommand=(vcmd, '%P'))
        self.ui_int.grid(column=1, row=2, padx=5, pady=5, sticky='nsew')
        self.ui_int.insert(0, '1.0')

        for i in range(2):
            self.ctl['spec_fr'].columnconfigure(i, weight=1)


    def add_CtlAutoFrame(self):
        """Draw the positioning controls."""

        self.cell_btn = {}

        for i,col in enumerate(CELL_NAMES):
            for j,row in enumerate(CELL_IDS):
                cell = col+row
                cell_id = CELL_NAMES[len(CELL_IDS)-j-1] + CELL_IDS[i]
                
                self.cell_btn[cell] = \
                    tk.Button(self.ctl['auto_fr'], text=cell, state='disabled',
                              command=lambda c=cell_id: self.ai.goto_cell(c))
                self.cell_btn[cell]\
                    .grid(column=i, row=len(CELL_IDS)-j-1,
                          sticky='nsew', padx=5, pady=5)
                
        ttk.Separator(self.ctl['auto_fr'], orient=tk.VERTICAL)\
           .grid(column=3, row=0, padx=5, pady=5, rowspan=3, sticky=(N,S))
        
        self.ctl['q_fr'] = tk.Frame(self.ctl['auto_fr'])
        ttk.Label(self.ctl['q_fr'], text='Number of chips in a cell')\
            .grid(column=0, row=0, padx=2, pady=5, sticky='e')
        self.ui_nchips = self.layout_CtlCombo('q_fr', 'NCHIPS',
                                              ['1', '2', '3'])
        self.ui_nchips.grid(column=1, row=0, padx=2, pady=5, sticky='w')

        self.ctl['q_fr'].columnconfigure(1, weight=1)

        self.ui_af = tk.Button(self.ctl['auto_fr'],
                               text='Auto-focus on sample',
                               state='disabled',
                               command=self.ai.auto_focus)
        self.ui_af.grid(column=4, row=0, sticky='nsew', padx=5, pady=5)
        
        self.ui_find = tk.Button(self.ctl['auto_fr'],
                                 text='Locate chip in a current cell',
                                 state='disabled',
                                 command=self.ai.find_chips)
        self.ui_find.grid(column=4, row=1, sticky='nsew', padx=5, pady=5)
        
        self.ctl['q_fr'].grid(column=4, row=2, sticky='nsew', padx=5, pady=5)
        self.ctl['auto_fr'].columnconfigure(3, weight=1)
        
            
    def add_CtlNaviFrame(self):
        """Draw the navigation controls."""

        lbls = ['+Y', '+Z', '-X', '+X', '-Y', '-Z']
        self.ui_mv = []
        
        for i,lbl in enumerate(lbls):
            self.ui_mv.append(tk.Button(self.ctl['navi_fr'], text=lbl,
                                        state='disabled',
                                        command=lambda l=lbl: self.ai.move(l)))
            self.ui_mv[-1]\
                .grid(column=i%2, row=i//2, sticky='nsew', padx=5, pady=5),

        ttk.Separator(self.ctl['navi_fr'], orient=tk.VERTICAL)\
           .grid(column=2, row=0, padx=5, pady=5, rowspan=3, sticky=(N,S))

        ttk.Label(self.ctl['navi_fr'], text='Horizontal step (mm)')\
           .grid(column=3, row=0, sticky='e', padx=5, pady=5)

        vcmd = self.ctl['navi_fr']\
                   .register(lambda s: self.val_num(s, 0.05, 10000))
        self.ui_hstep = ttk.Entry(self.ctl['navi_fr'],
                                  validate='focus',
                                  validatecommand=(vcmd, '%P'))
        self.ui_hstep.grid(column=4, row=0, sticky='w', padx=5, pady=5)
        self.ui_hstep.insert(0, '1.0')
        
        ttk.Label(self.ctl['navi_fr'], text='Vertical step (mm)')\
           .grid(column=3, row=1, sticky='e', padx=5, pady=5)

        self.ui_vstep = ttk.Entry(self.ctl['navi_fr'],
                                  validate='focus',
                                  validatecommand=(vcmd, '%P'))
        self.ui_vstep.grid(column=4, row=1, sticky='w', padx=5, pady=5)
        self.ui_vstep.insert(0, '1.0')

        self.gt_load = tk.Button(self.ctl['navi_fr'], text='Go to LOAD',
                                 state='disabled', command=self.ai.goto_load)
        self.gt_load.grid(column=3, row=2, sticky='nsew', padx=5, pady=5)
        self.gt_a1 = tk.Button(self.ctl['navi_fr'], text='Go to A1',
                               state='disabled', command=self.ai.goto_cell)
        self.gt_a1.grid(column=4, row=2, sticky='nsew', padx=5, pady=5)

        self.ctl['navi_fr'].columnconfigure(2, weight=1)


    def add_CtlCellSel(self):
        """Add the check button selection for which cells to scan."""
        
        self.cell_var, self.cell_chck = {}, {}

        for i,col in enumerate(CELL_NAMES):
            for j,row in enumerate(CELL_IDS):
                cell = col+row

                self.cell_var[cell] = tk.IntVar()
                self.cell_chck[cell] = \
                    tk.Checkbutton(self.ctl['cell_fr'], text=cell,
                                   variable=self.cell_var[cell])
                self.cell_chck[cell]\
                    .grid(column=i, row=len(CELL_IDS)-j-1,
                          sticky='nsew', padx=5, pady=5)

        for i in range(3):
            self.ctl['cell_fr'].columnconfigure(i, weight=1)

            
    def placeControls(self):
        """Create the frames for the control UI elements."""
        self.make_controlFrames()
        
        self.ui_conn = tk.Button(self.mf, text='CONNECT TO SPECTROMETER',
                                 command=self.connect, state=FEATURE_SET4)
        self.ui_conn.grid(column=0, row=0, columnspan=3, sticky='nsew',
                          padx=7, pady=10)

        self.add_CtlSettFrame()
        self.add_CtlSpecFrame()

        self.add_CtlNaviFrame()
        self.add_CtlAutoFrame()

        self.add_CtlCellSel()
        
        self.scan_btn = tk.Button(self.mf, text='SCAN CADDY',
                                  state='disabled',
                                  command=self.ai.scan_caddy)
        
        self.scan_btn.grid(column=1, row=9, columnspan=2, sticky='nsew',
                           padx=7, pady=10)

        tk.Button(self.mf, text='SAVE A1 POSITION', 
                  command=self.ai.reset_origin)\
          .grid(column=0, row=9, sticky='nsew', padx=7, pady=10)

        self.pos_lbl = tk.StringVar()
        ttk.Label(self.mf, text='')\
            .grid(column=0, row=10, columnspan=3, sticky='nsew')
        ttk.Label(self.mf, text='', textvariable=self.pos_lbl)\
            .grid(column=0, row=11, columnspan=3, sticky='nsew')

        self.mf.rowconfigure(10, weight=1)
        
        for i in range(3):
            self.mf.columnconfigure(i, weight=1)
            
