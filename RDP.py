import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import random
import math

DEFAULT_VALUES = {'A': 10000, 'E': 1000, 'I': 100, 'O': 10, 'U': 0, 'X': -10000}
REL_CHOICES = ['A', 'E', 'I', 'O', 'U', 'X']


class RDPApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title('Relationship Diagramming Process (RDP)')
        self.geometry('1200x700')

        self.values = DEFAULT_VALUES.copy()
        # Default to the provided 9x9 Relationship Chart (REL)
        self.n = 9
        self._last_n = self.n
        self.matrix = []

        self.sequence = []  # final Pi
        self.selection_steps = []  # logs for selection phase
        self.layout_steps = []

        self.current_phase = 'idle'  # 'selection' or 'layout'
        self.selection_index = 0
        self.layout_index = 0

        self.placed_positions = {}  # dept_id -> (x,y)
        self.started_for_matrix = False

        self._build_ui()

    def _build_ui(self):
        # Three column layout with user-resizable sashes
        paned = tk.PanedWindow(self, orient='horizontal')
        paned.pack(fill='both', expand=True, padx=8, pady=8)

        left = ttk.Frame(paned)
        middle = ttk.Frame(paned)
        right = ttk.Frame(paned)

        paned.add(left, minsize=180, stretch='never')
        paned.add(middle, minsize=360, stretch='always')
        paned.add(right, minsize=180, stretch='never')
        self._paned = paned
        self._in_paned_configure = False
        self._sashes_set = False
        self.after_idle(self._apply_initial_sashes)
        paned.bind('<Configure>', self._on_paned_configure)

        # Controls
        num_frame = tk.Frame(left, bg='#F1F4F8')
        num_frame.pack(fill='x', pady=(4, 8))
        tk.Label(num_frame, text='Number of Departments:', bg='#F1F4F8').pack(anchor='center', pady=(4, 2))
        self.n_var = tk.IntVar(value=self.n)
        n_spin = ttk.Spinbox(num_frame, from_=2, to=20, textvariable=self.n_var, width=6, command=self._on_n_change)
        n_spin.pack(anchor='center', pady=(0, 6))
        self.n_var.trace_add('write', lambda *args: self._on_n_change())

        rand_btn = tk.Button(left, text='Randomize Relations', command=self.randomize_matrix, bg='#ADD8E6')
        rand_btn.pack(fill='x', pady=4)

        cfg_btn = tk.Button(left, text='Configure Values', command=self.configure_values)
        cfg_btn.pack(fill='x', pady=4)

        start_btn = tk.Button(left, text='Start Calculation', command=self.start_calculation, bg='#90EE90')
        start_btn.pack(fill='x', pady=8)

        control_frame = ttk.Frame(left)
        control_frame.pack(fill='x', pady=8)
        # Row 1: Next / Jump / Run
        row1 = ttk.Frame(control_frame)
        row1.pack(fill='x')
        self.next_btn = tk.Button(row1, text='Next Step', command=self.next_step, state='disabled', bg='#00FFFF')
        self.next_btn.pack(side='left', expand=True, fill='x')
        self.jump_btn = tk.Button(row1, text='Jump to Layout Phase', command=self.jump_to_layout, bg='#FFA500', state='disabled')
        self.jump_btn.pack(side='left', expand=True, fill='x', padx=4)
        run_btn = tk.Button(row1, text='Run to End', command=self.run_to_end, state='disabled', bg='#90EE90')
        run_btn.pack(side='left', expand=True, fill='x')
        self.run_btn = run_btn
        # Row 2: Reset buttons
        row2 = ttk.Frame(control_frame)
        row2.pack(fill='x', pady=(6, 0))
        reset_btn = tk.Button(row2, text='Reset All', command=self.reset_all, bg='#FA8072')
        reset_btn.pack(side='left', expand=True, fill='x')
        reset_keep_btn = tk.Button(row2, text='Reset Steps', command=self.reset_steps_only, bg='#FAE28C')
        reset_keep_btn.pack(side='left', expand=True, fill='x', padx=4)

        preview_btn = tk.Button(left, text='Preview Report', command=self.preview_report, bg='#E0E0E0')
        preview_btn.pack(fill='x', pady=4)
        export_btn = tk.Button(left, text='Export Report (.txt)', command=self.export_report, bg='#D3D3D3')
        export_btn.pack(fill='x', pady=6)
        guide_btn = tk.Button(left, text='Guide', command=self.show_guide, bg='#E6E6FA')
        guide_btn.pack(fill='x', pady=4)
        info_btn = tk.Button(left, text='RDP Steps Info', command=self.show_steps_info, bg='#E6F3FF')
        info_btn.pack(fill='x', pady=4)

        ttk.Separator(left, orient='horizontal').pack(fill='x', pady=6)

        # Log area (left column, where matrix used to be)
        log_label = ttk.Label(left, text='Log:')
        log_label.pack(anchor='w')
        self.log_text = tk.Text(left, height=12)
        self.log_text.pack(fill='both', expand=True)

        # Middle: Matrix and Canvas (center column)
        self.matrix_container = ttk.Frame(middle)
        self.matrix_container.pack(fill='both', expand=True)
        self.matrix_canvas = tk.Canvas(self.matrix_container, height=360)
        self.matrix_scroll_y = ttk.Scrollbar(self.matrix_container, orient='vertical', command=self.matrix_canvas.yview)
        self.matrix_scroll_x = ttk.Scrollbar(self.matrix_container, orient='horizontal', command=self.matrix_canvas.xview)
        self.matrix_inner = ttk.Frame(self.matrix_canvas)
        self.matrix_inner_id = self.matrix_canvas.create_window((0,0), window=self.matrix_inner, anchor='nw')
        self.matrix_canvas.configure(yscrollcommand=self.matrix_scroll_y.set, xscrollcommand=self.matrix_scroll_x.set)
        self.matrix_canvas.pack(side='top', fill='both', expand=True)
        self.matrix_scroll_y.pack(side='right', fill='y')
        self.matrix_scroll_x.pack(side='bottom', fill='x')
        self.matrix_inner.bind('<Configure>', lambda e: self.matrix_canvas.configure(scrollregion=self.matrix_canvas.bbox('all')))
        self.matrix_canvas.bind('<Configure>', self._center_matrix_window)

        canvas_label = ttk.Label(middle, text='Layout Visualization:')
        canvas_label.pack(anchor='w')
        # Canvas with scrollbars
        canvas_frame = ttk.Frame(middle)
        canvas_frame.pack(fill='both', expand=True)
        self.canvas = tk.Canvas(canvas_frame, bg='white', height=280)
        self.canvas_h = ttk.Scrollbar(canvas_frame, orient='horizontal', command=self.canvas.xview)
        self.canvas_v = ttk.Scrollbar(canvas_frame, orient='vertical', command=self.canvas.yview)
        self.canvas.configure(xscrollcommand=self.canvas_h.set, yscrollcommand=self.canvas_v.set)
        self.canvas.grid(row=0, column=0, sticky='nsew')
        self.canvas_v.grid(row=0, column=1, sticky='ns')
        self.canvas_h.grid(row=1, column=0, sticky='ew')
        canvas_frame.rowconfigure(0, weight=1)
        canvas_frame.columnconfigure(0, weight=1)

        # storage for vars and widgets
        self.var_matrix = []
        self.om_matrix = []

        # Right column: Data tables (TCR top, WPV bottom)
        data_top = ttk.LabelFrame(right, text='TCR Table')
        data_top.pack(fill='both', expand=False)
        self.tcr_tree = ttk.Treeview(data_top, columns=('dept','tcr'), show='headings', height=8)
        self.tcr_tree.heading('dept', text='Dept')
        self.tcr_tree.heading('tcr', text='TCR')
        self.tcr_tree.pack(fill='both', expand=True)
        ttk.Label(data_top, text='Sequence (Pi):').pack(anchor='w', pady=(6, 0))
        self.sequence_text = tk.Text(data_top, height=2, wrap='word')
        self.sequence_text.pack(fill='x')
        self.sequence_text.config(state='disabled')

        data_bottom = ttk.LabelFrame(right, text='Placement Log & WPV')
        data_bottom.pack(fill='both', expand=True, pady=6)
        self.wpv_tree = ttk.Treeview(data_bottom, columns=('msg',), show='headings', height=8)
        self.wpv_tree.heading('msg', text='Placement Log & WPV')
        self.wpv_tree.pack(fill='both', expand=True)

        # initial matrix
        self.generate_matrix(init=True)

    def generate_matrix(self, init=True):
        try:
            n = int(self.n_var.get())
        except Exception:
            messagebox.showerror('Error', 'Invalid N')
            return
        if n < 2 or n > 20:
            messagebox.showerror('Error', 'N must be between 2 and 20')
            return
        self.n = n
        self._last_n = n
        self.started_for_matrix = False
        # initialize matrix with 'U' only if requested
        if init or not hasattr(self, 'matrix') or not self.matrix:
            self.matrix = [['U' for _ in range(n)] for _ in range(n)]
            for i in range(n):
                self.matrix[i][i] = 'U'

        # build widgets
        # reset any previous var/widget storage
        self.var_matrix = []
        self.om_matrix = []

        for child in self.matrix_inner.winfo_children():
            child.destroy()

        frame = ttk.Frame(self.matrix_inner)
        frame.pack()

        # header
        ttk.Label(frame, text='').grid(row=0, column=0)
        for j in range(n):
            ttk.Label(frame, text=f'D{j+1}').grid(row=0, column=j+1)
        for i in range(n):
            ttk.Label(frame, text=f'D{i+1}').grid(row=i+1, column=0)
            row_vars = []
            row_oms = []
            for j in range(n):
                display_val = '' if self.matrix[i][j] == 'U' else self.matrix[i][j]
                var = tk.StringVar(value=display_val)
                # use tk.OptionMenu so we can change background
                om = tk.OptionMenu(frame, var, *REL_CHOICES, command=lambda val, x=i, y=j: self._on_matrix_change(x, y, val))
                om.grid(row=i+1, column=j+1, padx=1, pady=1)
                if i == j:
                    om.config(state='disabled')
                row_vars.append(var)
                row_oms.append(om)
                # apply initial color
                self._apply_cell_color(i, j, om, display_val)
            self.var_matrix.append(row_vars)
            self.om_matrix.append(row_oms)

        self._clear_outputs()
        # enable jump button now that a matrix exists
        try:
            self.jump_btn.config(state='normal')
        except Exception:
            pass

    def _on_n_change(self):
        try:
            n = int(self.n_var.get())
        except Exception:
            return
        if n < 2 or n > 20:
            return
        if n == getattr(self, '_last_n', None):
            return
        self.generate_matrix(init=True)

    def _on_matrix_change(self, i, j, val):
        # val comes from OptionMenu; treat '' as 'U'
        chosen = val if val != '' else 'U'
        if 0 <= i < self.n and 0 <= j < self.n:
            self.started_for_matrix = False
            self._clear_outputs()
            # enforce symmetry
            self.matrix[i][j] = chosen
            self.matrix[j][i] = chosen
            # update displayed vars
            # set (i,j)
            disp = '' if chosen == 'U' else chosen
            try:
                self.var_matrix[i][j].set(disp)
                self._apply_cell_color(i, j, self.om_matrix[i][j], disp)
            except Exception:
                pass
            # set symmetric (j,i)
            try:
                self.var_matrix[j][i].set(disp)
                self._apply_cell_color(j, i, self.om_matrix[j][i], disp)
            except Exception:
                pass
            # calculations are only shown after Start Calculation

    def randomize_matrix(self):
        # Weighted randomization - sparser on strong relations
        weights = {
            'U': 50,
            'O': 20,
            'I': 15,
            'E': 10,
            'A': 5,
            'X': 0.5
        }
        letters = list(weights.keys())
        w = [weights[l] for l in letters]
        for i in range(self.n):
            for j in range(i, self.n):
                if i == j:
                    self.matrix[i][j] = 'U'
                else:
                    val = random.choices(letters, weights=w, k=1)[0]
                    # Very small chance for X, ensure less than 5% overall
                    self.matrix[i][j] = val
                    self.matrix[j][i] = val
        # rebuild widgets from existing self.matrix without re-initializing it
        self.generate_matrix(init=False)

    def configure_values(self):
        dlg = ConfigValuesDialog(self, self.values)
        self.wait_window(dlg)
        if dlg.result:
            self.values = dlg.result
            self.started_for_matrix = False
            self.log(f'Values updated: {self.values}')
            # refresh TCR view if present
            self.update_tcr_view()

    def start_calculation(self):
        if self.started_for_matrix:
            return
        # compute selection sequence following the strict logic
        self.reset_state_for_run()
        self.compute_tcr()
        self.build_selection_sequence()
        self.current_phase = 'selection'
        self.selection_index = 0
        self.log('Selection phase prepared. Use Next Step, Jump to Layout, or Run to End.')
        self.log('Next Step: advances one step at a time.')
        self.log('Jump to Layout: completes selection and prepares layout.')
        self.log('Run to End: completes all steps quickly.')
        self.next_btn.config(state='normal')
        self.run_btn.config(state='normal')
        try:
            self.jump_btn.config(state='normal')
        except Exception:
            pass
        self.started_for_matrix = True

    def jump_to_layout(self):
        # Execute Steps 1-6 (calculation & sequencing) and prepare layout, stop before placements
        # run silently and enable layout controls
        self.reset_state_for_run()
        self.compute_tcr()
        self.build_selection_sequence()
        for s in self.selection_steps:
            self.log(s)
        # prepare layout steps (but do not execute any placement)
        self.prepare_layout()
        # move to layout phase, but don't show ghosts or place anything
        self.current_phase = 'layout'
        self.layout_index = 0
        # Log summary: show final TCR and full sequence
        self.log('Jumped to Layout Phase: Steps 1-6 completed (stopped before first placement).')
        self.log('Final Sequence: ' + ' -> '.join(f'D{d+1}' for d in self.sequence))
        # populate TCR view already done by compute_tcr; ensure WPV is cleared
        try:
            self.wpv_tree.delete(*self.wpv_tree.get_children())
        except Exception:
            pass
        # enable flow controls for layout continuation
        self.next_btn.config(state='normal')
        self.run_btn.config(state='normal')
        # disable jump to avoid repeated jumps
        try:
            self.jump_btn.config(state='disabled')
        except Exception:
            pass

    def reset_state_for_run(self):
        self.sequence = []
        self.selection_steps = []
        self.layout_steps = []
        self.placed_positions = {}
        self.canvas.delete('all')
        self.selection_index = 0
        self.layout_index = 0
        self.layout_showing_ghosts = False

    def compute_tcr(self, log_result=False):
        # TCR excludes 'X' relationships
        self.tcr = [0 for _ in range(self.n)]
        for i in range(self.n):
            s = 0
            for j in range(self.n):
                if i == j:
                    continue
                rel = self.matrix[i][j]
                if rel == 'X':
                    continue
                s += self.values.get(rel, 0)
            self.tcr[i] = s
        if log_result:
            self.log(f'Computed TCR: {self.tcr}')
        self.update_tcr_view()

    def update_tcr_view(self):
        # populate treeview
        try:
            self.tcr_tree.delete(*self.tcr_tree.get_children())
        except Exception:
            return
        for i in range(self.n):
            t = self.tcr[i] if hasattr(self, 'tcr') else 0
            self.tcr_tree.insert('', 'end', values=(f'D{i+1}', t))

    def _apply_cell_color(self, i, j, widget, display_val):
        # display_val is '' for U or letter
        letter = display_val if display_val != '' else 'U'
        color_map = {
            'A': '#006400',
            'E': '#228B22',
            'I': '#90EE90',
            'O': '#FFFACD',
            'U': '#FFFFFF',
            'X': '#FF4C4C'
        }
        bg = color_map.get(letter, '#FFFFFF')
        # choose text color
        text_col = '#ffffff' if letter in ('A', 'E', 'X') else '#000000'
        try:
            widget.config(bg=bg, fg=text_col)
            # menu colors
            menu = widget['menu']
            menu.config(bg=bg, fg=text_col)
        except Exception:
            try:
                widget.configure(background=bg)
            except Exception:
                pass

    def build_selection_sequence(self):
        P = set(range(self.n))
        Pi = [None] * self.n
        idx_front = 0
        idx_back = self.n - 1

        # Step 4: select first by highest TCR
        while P:
            if idx_front == 0:
                # select first
                max_tcr = max((self.tcr[i] for i in P))
                candidates = [i for i in P if self.tcr[i] == max_tcr]
                if len(candidates) > 1:
                    # tie-breaker by count of 'A'
                    best = []
                    best_count = -1
                    for c in candidates:
                        countA = sum(1 for j in range(self.n) if self.matrix[c][j] == 'A')
                        if countA > best_count:
                            best = [c]
                            best_count = countA
                        elif countA == best_count:
                            best.append(c)
                    if len(best) > 1:
                        chosen = random.choice(best)
                    else:
                        chosen = best[0]
                else:
                    chosen = candidates[0]
                Pi[idx_front] = chosen
                P.remove(chosen)
                self.selection_steps.append(f'Step4: Dept D{chosen+1} selected (TCR={self.tcr[chosen]}). Placed at front index {idx_front}.')
                idx_front += 1
                current = chosen
                # handle X relationships
                while True:
                    x_related = [p for p in list(P) if self.matrix[p][current] == 'X']
                    if not x_related:
                        break
                    for xr in x_related:
                        Pi[idx_back] = xr
                        P.remove(xr)
                        self.selection_steps.append(f'Step5: Dept D{xr+1} has X with D{current+1}; placed at back index {idx_back}.')
                        idx_back -= 1
                continue

            # now Step 6: select next based on highest relation to assigned front part
            assigned = [x for x in Pi[:idx_front] if x is not None]
            if not assigned:
                break
            # priority order list
            priority = ['A', 'E', 'I', 'O']
            candidate_priority = {}
            for p in P:
                best_rank = None
                for r in priority:
                    # if p has relation r with any of assigned
                    if any(self.matrix[p][a] == r or self.matrix[a][p] == r for a in assigned):
                        best_rank = r
                        break
                candidate_priority[p] = best_rank

            # find highest priority present
            chosen_candidates = [p for p, r in candidate_priority.items() if r is not None]
            if chosen_candidates:
                # group by rank according to priority order
                for r in priority:
                    group = [p for p in chosen_candidates if candidate_priority[p] == r]
                    if group:
                        # tie-breaker by highest TCR
                        best_tcr = max(self.tcr[g] for g in group)
                        bests = [g for g in group if self.tcr[g] == best_tcr]
                        chosen = random.choice(bests) if len(bests) > 1 else bests[0]
                        Pi[idx_front] = chosen
                        P.remove(chosen)
                        self.selection_steps.append(f'Step6: Dept D{chosen+1} selected with relation {r} to assigned set; placed at front index {idx_front}.')
                        idx_front += 1
                        current = chosen
                        # handle X related to this current
                        while True:
                            x_related = [p for p in list(P) if self.matrix[p][current] == 'X']
                            if not x_related:
                                break
                            for xr in x_related:
                                Pi[idx_back] = xr
                                P.remove(xr)
                                self.selection_steps.append(f'Step5: Dept D{xr+1} has X with D{current+1}; placed at back index {idx_back}.')
                                idx_back -= 1
                        break
                continue

            # if no one has any A/E/I/O with assigned, simply pick highest TCR remaining
            max_tcr = max((self.tcr[i] for i in P))
            candidates = [i for i in P if self.tcr[i] == max_tcr]
            chosen = random.choice(candidates)
            Pi[idx_front] = chosen
            P.remove(chosen)
            self.selection_steps.append(f'Fallback: Dept D{chosen+1} selected by TCR; placed at front index {idx_front}.')
            idx_front += 1

        # fill any remaining None (shouldn't happen)
        for i in range(self.n):
            if Pi[i] is None:
                remaining = [x for x in range(self.n) if x not in Pi]
                if remaining:
                    Pi[i] = remaining.pop(0)

        self.sequence = Pi
        self._update_sequence_view()

    def next_step(self):
        if self.current_phase == 'selection':
            if self.selection_index < len(self.selection_steps):
                s = self.selection_steps[self.selection_index]
                self.log(s)
                self.selection_index += 1
                # if finished selection steps, switch to layout
                if self.selection_index >= len(self.selection_steps):
                    self.current_phase = 'layout'
                    self.prepare_layout()
                    self.log('Switching to layout phase. Use Next Step to place departments.')
                return
        if self.current_phase == 'layout':
            if self.layout_index < len(self.layout_steps):
                step = self.layout_steps[self.layout_index]
                # unpack (dept, best_spot, best_wpv, details, candidates)
                if len(step) == 5:
                    dept, best_spot, wpv, details, candidates = step
                else:
                    dept, best_spot, wpv, details = step
                    candidates = []

                # If ghosts not yet shown for this step, draw them and wait for confirmation
                if not getattr(self, 'layout_showing_ghosts', False):
                    # draw candidate ghost blocks with WPV overlay
                    self._draw_ghosts(candidates)
                    # populate WPV tree with candidates for debug
                    try:
                        self.wpv_tree.delete(*self.wpv_tree.get_children())
                    except Exception:
                        pass
                    for spot, wpv_val, det in candidates:
                        self.wpv_tree.insert('', 'end', values=(f'D{dept+1} candidate {spot}: WPV={wpv_val:.1f}',))
                    self.log(f'Evaluating placement for D{dept+1}: showing {len(candidates)} candidate spots with WPV')
                    self.layout_showing_ghosts = True
                    return

                # confirmation: remove ghosts and place the department
                self.canvas.delete('ghost')
                self.placed_positions[dept] = best_spot
                self.draw_layout()
                self.log(details)
                # record placement in WPV log
                try:
                    self.wpv_tree.insert('', 'end', values=(f'D{dept+1}: Placed at {best_spot} with WPV={wpv:.1f}',))
                except Exception:
                    pass
                self.layout_index += 1
                self.layout_showing_ghosts = False
                if self.layout_index >= len(self.layout_steps):
                    self.log('Layout complete.')
                    self.next_btn.config(state='disabled')
                    self.run_btn.config(state='disabled')
                return

    def run_to_end(self):
        if self.current_phase == 'selection':
            # output remaining selection logs
            for i in range(self.selection_index, len(self.selection_steps)):
                self.log(self.selection_steps[i])
            self.selection_index = len(self.selection_steps)
            self.current_phase = 'layout'
            self.prepare_layout()
            self.log('Switching to layout phase.')

        if self.current_phase == 'layout':
            # perform all layout placements
            for i in range(self.layout_index, len(self.layout_steps)):
                step = self.layout_steps[i]
                if len(step) == 5:
                    dept, coord, wpv, details, candidates = step
                else:
                    dept, coord, wpv, details = step
                self.placed_positions[dept] = coord
                self.draw_layout()
                self.log(details)
                # add to WPV log
                try:
                    self.wpv_tree.insert('', 'end', values=(f'D{dept+1}: Placed at {coord} with WPV={wpv:.1f}',))
                except Exception:
                    pass
            self.layout_index = len(self.layout_steps)
            self.log('Layout complete.')
            self.next_btn.config(state='disabled')
            self.run_btn.config(state='disabled')

    def prepare_layout(self):
        # build layout steps from self.sequence
        # Each layout_steps entry: (dept, best_spot, best_wpv, best_details, candidates)
        self.layout_steps = []
        center = (0, 0)
        occupied = {}
        for k, dept in enumerate(self.sequence):
            if k == 0:
                coord = center
                wpv = 0
                details = f'Place D{dept+1} at center {coord} (seed).'
                self.layout_steps.append((dept, coord, wpv, details, []))
                occupied[coord] = dept
                continue

            # compute available spots adjacent to cluster
            available = set()
            xs = [p[0] for p in occupied.keys()]
            ys = [p[1] for p in occupied.keys()]
            center_x = sum(xs)/len(xs)
            center_y = sum(ys)/len(ys)
            for (ox, oy) in occupied.keys():
                for dx in [-1, 0, 1]:
                    for dy in [-1, 0, 1]:
                        if dx == 0 and dy == 0:
                            continue
                        cand = (ox+dx, oy+dy)
                        if cand in occupied:
                            continue
                        available.add(cand)

            if not available:
                # expand around center
                available.add((0, -1))

            # sort available starting from Western Edge and moving CCW
            def ang_key(pt):
                ax = pt[0] - center_x
                ay = pt[1] - center_y
                angle = math.atan2(ay, ax)
                # shift so that west (pi or -pi) -> 0 and increase CCW
                key = (angle - math.pi) % (2*math.pi)
                # also use distance as secondary key (closer first)
                dist = math.hypot(ax, ay)
                return (key, dist)

            sorted_avail = sorted(available, key=ang_key)

            # evaluate WPV for each spot
            best_spot = None
            best_wpv = -1e18
            best_details = ''
            candidates = []
            for spot in sorted_avail:
                wpv = 0.0
                detail_parts = []
                # store candidate wpv for debug
                # neighbor contributions
                sx, sy = spot
                # look at neighbors of this spot that are occupied
                for ox, oy in occupied.keys():
                    nx = ox - sx
                    ny = oy - sy
                    # determine adjacency
                    if abs(nx) + abs(ny) == 1 and (abs(nx)+abs(ny)) != 0:
                        weight = 1.0
                    elif abs(nx) == 1 and abs(ny) == 1:
                        weight = 0.5
                    else:
                        continue
                    neighbor_dept = occupied[(ox, oy)]
                    rel_letter = self.matrix[dept][neighbor_dept]
                    rel_value = self.values.get(rel_letter, 0)
                    contrib = rel_value * weight
                    wpv += contrib
                    detail_parts.append(f'with D{neighbor_dept+1}({rel_letter})*{weight}={contrib}')
                details = f'Place D{dept+1} at {spot}: WPV={wpv}. ' + '; '.join(detail_parts)
                if wpv > best_wpv:
                    best_wpv = wpv
                    best_spot = spot
                    best_details = details
                # append candidate
                candidates.append((spot, wpv, details))
            # tie-breaker already handled by sorted_avail order because we iterate in that order and only replace on strictly greater
            if best_spot is None:
                best_spot = sorted_avail[0]
                best_details = f'Place D{dept+1} at {best_spot} by default.'
            # record and mark occupied
            self.layout_steps.append((dept, best_spot, best_wpv, best_details, candidates))
            occupied[best_spot] = dept

    def draw_layout(self):
        # only remove placed items; preserve ghosts
        self.canvas.delete('placed')
        if not self.placed_positions:
            return
        cell = 60
        # canvas center coordinates
        cw = max(self.canvas.winfo_width(), 200)
        ch = max(self.canvas.winfo_height(), 200)
        center_px = cw/2
        center_py = ch/2

        # compute bounding extents to set scrollregion
        pxs = []
        pys = []
        for dept, (x,y) in self.placed_positions.items():
            px = center_px + x*cell
            py = center_py - y*cell
            pxs.append(px)
            pys.append(py)

        min_px = min(pxs) - cell
        max_px = max(pxs) + cell
        min_py = min(pys) - cell
        max_py = max(pys) + cell
        self.canvas.config(scrollregion=(min_px, min_py, max_px, max_py))

        # draw each placed
        for dept, (x,y) in self.placed_positions.items():
            px = center_px + x*cell
            py = center_py - y*cell
            x1 = px - cell/2 + 5
            y1 = py - cell/2 + 5
            x2 = px + cell/2 - 5
            y2 = py + cell/2 - 5
            color = f'#{(hash(dept)%8+1)*25:02x}{(hash(dept+3)%8+2)*20:02x}{(hash(dept+5)%8+3)*18:02x}'
            self.canvas.create_rectangle(x1,y1,x2,y2, fill=color, tags=('placed',))
            self.canvas.create_text((x1+x2)/2, (y1+y2)/2, text=f'D{dept+1}', tags=('placed',))

    def grid_to_canvas(self, x, y, cell=60):
        cw = max(self.canvas.winfo_width(), 200)
        ch = max(self.canvas.winfo_height(), 200)
        center_px = cw/2
        center_py = ch/2
        px = center_px + x*cell
        py = center_py - y*cell
        x1 = px - cell/2 + 5
        y1 = py - cell/2 + 5
        x2 = px + cell/2 - 5
        y2 = py + cell/2 - 5
        return x1, y1, x2, y2, px, py

    def _draw_ghosts(self, candidates):
        # candidates: list of (spot, wpv, details)
        self.canvas.delete('ghost')
        for spot, wpv, details in candidates:
            x1, y1, x2, y2, px, py = self.grid_to_canvas(spot[0], spot[1])
            self.canvas.create_rectangle(x1, y1, x2, y2, outline='#555555', dash=(3,3), tags=('ghost',))
            self.canvas.create_text(px, py, text=f'{wpv:.1f}', fill='#000000', tags=('ghost',))

    def reset_all(self):
        # Reset algorithmic state
        self.reset_state_for_run()
        self.started_for_matrix = False
        self.sequence = []
        self.selection_steps = []
        self.layout_steps = []
        self.current_phase = 'idle'
        self.selection_index = 0
        self.layout_index = 0
        self.placed_positions = {}

        self._clear_outputs()

        # Reset matrix inputs to default 'U' and update widgets
        try:
            if hasattr(self, 'matrix') and self.matrix:
                for i in range(self.n):
                    for j in range(self.n):
                        self.matrix[i][j] = 'U'
                        try:
                            # set display var to '' for U
                            if self.var_matrix and self.var_matrix[i][j] is not None:
                                self.var_matrix[i][j].set('')
                            if self.om_matrix and self.om_matrix[i][j] is not None:
                                self._apply_cell_color(i, j, self.om_matrix[i][j], '')
                        except Exception:
                            pass
        except Exception:
            pass

        # disable flow controls until next calculation
        try:
            self.next_btn.config(state='disabled')
            self.run_btn.config(state='disabled')
        except Exception:
            pass
        # re-enable jump (matrix may still be present) but default to disabled until generation
        try:
            self.jump_btn.config(state='disabled')
        except Exception:
            pass

    def reset_steps_only(self):
        # Keep matrix; clear progression, logs, and placements
        self.reset_state_for_run()
        self.current_phase = 'idle'
        self.selection_index = 0
        self.layout_index = 0
        self.started_for_matrix = False
        self._clear_outputs()
        try:
            self.next_btn.config(state='disabled')
            self.run_btn.config(state='disabled')
            self.jump_btn.config(state='normal')
        except Exception:
            pass

    def log(self, msg):
        self.log_text.insert('end', f'- {msg}\n')
        self.log_text.see('end')

    def _clear_outputs(self):
        try:
            self.canvas.delete('all')
        except Exception:
            pass
        try:
            self.log_text.delete('1.0', 'end')
        except Exception:
            pass
        try:
            self.tcr_tree.delete(*self.tcr_tree.get_children())
        except Exception:
            pass
        try:
            self.wpv_tree.delete(*self.wpv_tree.get_children())
        except Exception:
            pass
        try:
            self._update_sequence_view(clear=True)
        except Exception:
            pass

    def show_guide(self):
        guide = (
            'Quick Guide\n'
            '- Set the number of departments; the matrix updates automatically.\n'
            '- Edit relationships or use Randomize Relations.\n'
            '- Start Calculation: computes TCR and prepares selection steps.\n'
            '- Next Step: advances one step at a time.\n'
            '- Jump to Layout: completes selection and prepares layout.\n'
            '- Run to End: completes all steps quickly.\n'
            '- Reset Steps: keeps the matrix, clears progress.\n'
            '- Reset All: clears everything.\n'
            '- Export Report: exports inputs and the algorithm flow.'
        )
        messagebox.showinfo('Guide', guide)

    def show_steps_info(self):
        info = (
            'RDP Algorithm Steps (Summary)\n'
            '- Step 1: define departments and relationship values (A/E/I/O/U/X).\n'
            '- Step 2: build the relationship matrix (REL).\n'
            '- Step 3: compute TCR for each department and initialize the sequence.\n'
            '- Step 4: pick the first department with the highest TCR.\n'
            '- Step 5: if a department has an X relation to the current, place it at the back.\n'
            '- Step 6: choose the next department with the highest-priority relation\n'
            '  to the already selected set (A > E > I > O).\n'
            '- Fallback: if no A/E/I/O exists, pick by highest TCR.\n'
            '\n'
            'Note: The log shows which rule was applied, not a global step count.'
        )
        messagebox.showinfo('RDP Steps Info', info)

    def _update_sequence_view(self, clear=False):
        try:
            self.sequence_text.config(state='normal')
            self.sequence_text.delete('1.0', 'end')
            if clear or not self.sequence:
                self.sequence_text.insert('end', '-')
            else:
                seq = ' -> '.join(f'D{d+1}' for d in self.sequence)
                self.sequence_text.insert('end', seq)
            self.sequence_text.config(state='disabled')
        except Exception:
            pass

    def _apply_initial_sashes(self):
        if not self._sashes_set:
            self._set_paned_sashes()

    def _on_paned_configure(self, event):
        if self._in_paned_configure:
            return
        self._set_paned_sashes()

    def _set_paned_sashes(self):
        try:
            self._in_paned_configure = True
            total = self._paned.winfo_width()
            if total <= 0:
                return
            left_w = max(180, int(total * 0.20))
            right_start = max(left_w + 360, int(total * 0.80))
            self._paned.sash_place(0, left_w, 0)
            self._paned.sash_place(1, right_start, 0)
            self._sashes_set = True
        finally:
            self._in_paned_configure = False

    def _center_matrix_window(self, event):
        try:
            self.matrix_canvas.update_idletasks()
            cw = self.matrix_canvas.winfo_width()
            ch = self.matrix_canvas.winfo_height()
            iw = self.matrix_inner.winfo_reqwidth()
            ih = self.matrix_inner.winfo_reqheight()
            x = max((cw - iw) // 2, 0)
            y = max((ch - ih) // 2, 0)
            self.matrix_canvas.coords(self.matrix_inner_id, x, y)
        except Exception:
            pass

    def _compute_tcr_values(self):
        # compute TCR without mutating UI/logs
        tcr = [0 for _ in range(self.n)]
        for i in range(self.n):
            s = 0
            for j in range(self.n):
                if i == j:
                    continue
                rel = self.matrix[i][j]
                if rel == 'X':
                    continue
                s += self.values.get(rel, 0)
            tcr[i] = s
        return tcr

    def _layout_ascii(self):
        if not self.placed_positions:
            return 'No layout placements yet.'
        xs = [p[0] for p in self.placed_positions.values()]
        ys = [p[1] for p in self.placed_positions.values()]
        min_x, max_x = min(xs), max(xs)
        min_y, max_y = min(ys), max(ys)
        cols = max_x - min_x + 1
        rows = max_y - min_y + 1
        cell_w = 4
        grid = [[' ' * cell_w for _ in range(cols)] for _ in range(rows)]
        for dept, (x, y) in self.placed_positions.items():
            gx = x - min_x
            gy = max_y - y
            label = f'D{dept+1}'
            label = label[:cell_w].ljust(cell_w)
            if 0 <= gy < rows and 0 <= gx < cols:
                grid[gy][gx] = label
        lines = [''.join(row).rstrip() for row in grid]
        header = f'Layout grid (x={min_x}..{max_x}, y={min_y}..{max_y}):'
        return header + '\n' + '\n'.join(lines)

    def _build_report_lines(self):
        lines = []
        lines.append('RDP Report')
        lines.append('=' * 40)
        lines.append('1) Inputs')
        lines.append(f'- Number of Departments: {self.n}')
        lines.append('- Relationship Values:')
        for k in REL_CHOICES:
            lines.append(f'  - {k} = {self.values.get(k, 0)}')
        lines.append('')
        lines.append('2) Relationship Matrix')
        header = '     ' + ' '.join([f'D{j+1:02d}' for j in range(self.n)])
        lines.append(header)
        for i in range(self.n):
            row = [self.matrix[i][j] if self.matrix[i][j] != '' else 'U' for j in range(self.n)]
            lines.append(f'D{i+1:02d}  ' + '  '.join(row))
        lines.append('')
        lines.append('3) TCR')
        if self.started_for_matrix and hasattr(self, 'tcr'):
            tcr_vals = self._compute_tcr_values()
            lines.append('- ' + ', '.join([f'D{i+1}={tcr_vals[i]}' for i in range(self.n)]))
        else:
            lines.append('- (Not computed yet)')
        lines.append('')
        lines.append('4) Final Sequence')
        if self.started_for_matrix and self.sequence:
            lines.append('- ' + ' -> '.join([f'D{d+1}' for d in self.sequence]))
        else:
            lines.append('- (Not computed yet)')
        lines.append('')
        lines.append('5) Selection Steps')
        if self.started_for_matrix and self.selection_steps:
            for s in self.selection_steps:
                if s.startswith('Step4:'):
                    lines.append('- Step 4: ' + s.replace('Step4: ', ''))
                elif s.startswith('Step5:'):
                    lines.append('- Step 5: ' + s.replace('Step5: ', ''))
                elif s.startswith('Step6:'):
                    lines.append('- Step 6: ' + s.replace('Step6: ', ''))
                elif s.startswith('Fallback:'):
                    lines.append('- Fallback: ' + s.replace('Fallback: ', ''))
                else:
                    lines.append('- ' + s)
        else:
            lines.append('- (Not computed yet)')
        lines.append('')
        lines.append('6) Layout Steps')
        if self.started_for_matrix and self.layout_steps:
            for step in self.layout_steps:
                if len(step) >= 4:
                    dept, best_spot, wpv, details = step[:4]
                    lines.append(f'- D{dept+1} at {best_spot} WPV={wpv:.1f}')
                    lines.append(f'  - {details}')
        else:
            lines.append('- (Not computed yet)')
        lines.append('')
        lines.append('7) Final Layout')
        if self.placed_positions:
            lines.append('- Positions:')
            for dept in sorted(self.placed_positions.keys()):
                lines.append(f'  - D{dept+1}: {self.placed_positions[dept]}')
        else:
            lines.append('- (No placements yet)')
        lines.append(self._layout_ascii())
        lines.append('')
        # omit UI logs; report focuses on algorithm flow only
        return lines

    def preview_report(self):
        lines = self._build_report_lines()
        preview = '\n'.join(lines)
        win = tk.Toplevel(self)
        win.title('Report Preview')
        win.geometry('800x600')
        win.transient(self)

        frame = ttk.Frame(win)
        frame.pack(fill='both', expand=True, padx=8, pady=8)

        text = tk.Text(frame, wrap='none')
        y_scroll = ttk.Scrollbar(frame, orient='vertical', command=text.yview)
        x_scroll = ttk.Scrollbar(frame, orient='horizontal', command=text.xview)
        text.configure(yscrollcommand=y_scroll.set, xscrollcommand=x_scroll.set)

        text.grid(row=0, column=0, sticky='nsew')
        y_scroll.grid(row=0, column=1, sticky='ns')
        x_scroll.grid(row=1, column=0, sticky='ew')
        frame.rowconfigure(0, weight=1)
        frame.columnconfigure(0, weight=1)

        text.insert('1.0', preview)
        text.config(state='disabled')

    def export_report(self):
        path = filedialog.asksaveasfilename(
            title='Export RDP Report',
            defaultextension='.txt',
            filetypes=[('Text files', '*.txt'), ('All files', '*.*')]
        )
        if not path:
            return

        lines = self._build_report_lines()
        try:
            with open(path, 'w', encoding='utf-8') as f:
                f.write('\n'.join(lines))
            messagebox.showinfo('Export', f'Report saved to:\n{path}')
        except Exception as e:
            messagebox.showerror('Export Error', f'Failed to save report:\n{e}')


class ConfigValuesDialog(tk.Toplevel):
    def __init__(self, parent, values):
        super().__init__(parent)
        self.title('Configure Relationship Values')
        self.result = None
        self.transient(parent)
        self.grab_set()

        self.vars = {}
        row = 0
        for key in REL_CHOICES:
            ttk.Label(self, text=key).grid(row=row, column=0, padx=6, pady=6)
            var = tk.StringVar(value=str(values.get(key, 0)))
            self.vars[key] = var
            ttk.Entry(self, textvariable=var).grid(row=row, column=1, padx=6, pady=6)
            row += 1
        btn = ttk.Button(self, text='OK', command=self.on_ok)
        btn.grid(row=row, column=0, columnspan=2, pady=8)

    def on_ok(self):
        res = {}
        try:
            for k, v in self.vars.items():
                res[k] = float(v.get()) if '.' in v.get() else int(v.get())
        except Exception:
            messagebox.showerror('Error', 'Invalid numeric value')
            return
        self.result = res
        self.destroy()


if __name__ == '__main__':
    app = RDPApp()
    app.mainloop()
