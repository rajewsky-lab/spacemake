import base64
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import logging
import scanpy as sc
import functools
import uuid
import math
import anndata as ad

from itertools import cycle
from scipy.stats import gaussian_kde
from dataclasses import dataclass, field
from typing import Dict, List, Any, Optional, Callable
from io import BytesIO

cpalette = {
    "grey": "#999999",
    "light_orange": "#E69F00",
    "light_blue": "#56B4E9",
    "green": "#009E73",
    "yellow": "#F0E442",
    "blue": "#0072B2",
    "orange": "#D55E00",
    "pink": "#CC79A7",
}

metrics_colors = {
    "umis": "#E69F00",
    "genes": "#56B4E9",
    "reads": "#009E73",
    "pcr": "#CC79A7",
    "pct_counts_mt": "black",
}

nucl_clrs = {
    "A": "#F5C900",
    "C": "#F55D59",
    "T": "#3AA861",
    "G": "#7772F5",
    "N": "#999999",
}

parula_dict = {
    1: "#352a87",
    2: "#2058b0",
    3: "#1f7eb8",
    4: "#28a7d7",
    5: "#38d7e3",
    6: "#99d4d0",
    7: "#aacca1",
    8: "#bbcc74",
    9: "#cbcc49",
    10: "#e0d317",
}

PCT_DOWNSAMPLE_TO_PLOT = [20, 40, 60, 80, 100]

logger_name = "spacemake.pl"
logger = logging.getLogger(logger_name)

def ensure_anndata(func):
    """
    Decorator to ensure the first argument is an AnnData object.
    If a string is provided, it will load the file as AnnData.
    """
    @functools.wraps(func)
    def wrapper(adata, *args, **kwargs):
        if isinstance(adata, str):
            adata = ad.read_h5ad(adata)
        return func(adata, *args, **kwargs)
    return wrapper

@dataclass
class Plot:
    """Represents a single plot with its metadata."""
    title: str
    description: str
    plot_func: Callable
    
    def generate(self) -> str:
        """Generate the HTML for this plot."""
        return f"""
        <div class="plot-section">
            <h4>{self.title}</h4>
            <p class="text-muted">{self.description}</p>
            {self._get_plot_html()}
        </div>
        """
    
    def _get_plot_html(self) -> str:
        """Convert matplotlib figure to HTML."""
        try:
            with plt.ioff():
                result = self.plot_func()
                if result is None:
                    return '<div class="alert alert-warning">No plot data available</div>'
                
                if isinstance(result, tuple) and len(result) == 2:
                    fig, ax = result
                else:
                    fig = result
                
                buf = BytesIO()
                fig.savefig(buf, format='png', bbox_inches='tight', dpi=150, facecolor='white')
                plt.close(fig)
                
                data = base64.b64encode(buf.getvalue()).decode('utf-8')
                return f'<div class="text-center"><img src="data:image/png;base64,{data}" class="img-fluid" style="max-width: 100%; height: auto;"/></div>'
        except Exception as e:
            # just empty if plot cannot be shown
            return f'<div class="alert alert-primary">Plot not available for this configuration</div>'
    
@dataclass
class PlotGroup:
    """Represents a group of related plots."""
    name: str
    description: str
    plots: List[Plot]
    
    def get_clean_id(self) -> str:
        """Get a clean ID for HTML use."""
        return "".join(c.lower() if c.isalnum() else "-" for c in self.name).strip("-")
    
    def generate(self, active: bool = False) -> str:
        """Generate HTML for all plots in the group."""
        plots_html = "\n".join(plot.generate() for plot in self.plots)
        clean_id = self.get_clean_id()
        
        return f"""
        <div class="tab-pane fade {'show active' if active else ''}" 
             id="{clean_id}" 
             role="tabpanel" 
             aria-labelledby="{clean_id}-tab">
            <h3>{self.name}</h3>
            <p class="lead text-muted">{self.description}</p>
            {plots_html}
        </div>
        """

@dataclass
class Column:
    """Represents a single column in a table."""
    name: str
    description: str
    formatter: Optional[Callable[[Any], str]] = None
    
    def format_value(self, value: Any) -> str:
        """Format a value according to the column's rules."""
        if self.formatter is not None:
            return self.formatter(value)
        return str(value)

@dataclass
class TableStyle:
    """Defines the styling for a table."""
    table_classes: List[str] = field(default_factory=lambda: ["table", "table-striped"])
    container_classes: List[str] = field(default_factory=lambda: ["table-container"])
    header_classes: List[str] = field(default_factory=lambda: ["thead-light"])
    row_classes: List[str] = field(default_factory=lambda: [])
    
    def get_table_class(self) -> str:
        return " ".join(self.table_classes)
    
    def get_container_class(self) -> str:
        return " ".join(self.container_classes)
    
    def get_header_class(self) -> str:
        return " ".join(self.header_classes)
    
    def get_row_class(self) -> str:
        return " ".join(self.row_classes)

@dataclass
class DataFrameTable:
    """A table based on a pandas DataFrame with custom column formatting."""
    data: pd.DataFrame
    title: str
    description: str
    style: TableStyle = field(default_factory=TableStyle)
    columns: Optional[Dict[str, Column]] = None
    
    def __post_init__(self):
        """Initialize columns if not provided."""
        if self.columns is None:
            self.columns = {
                col: Column(
                    name=str(col),
                    description=str(col),
                    formatter=None
                ) for col in self.data.columns
            }
    
    def generate(self):
        from IPython.display import HTML

        table_html = f"""
        <div class="table-section">
            <h4>{self.title}</h4>
            <p class="text-muted">{self.description}</p>
            <div class="{self.style.get_container_class()}">
                <table class="{self.style.get_table_class()}">
                    <thead class="{self.style.get_header_class()}">
                        <tr>
        """
        
        # Add headers
        for col_name, col in self.columns.items():
            table_html += f'<th scope="col" title="{col.description}">{col.name}</th>'
        
        table_html += """
                        </tr>
                    </thead>
                    <tbody>
        """

        display_data = self.data.reset_index()
        
        # Add data rows
        for _, row in display_data.iterrows():
            table_html += f'<tr class="{self.style.get_row_class()}">'
            for col_name, col in self.columns.items():
                formatted_value = col.format_value(row[col_name])
                table_html += f"<td>{formatted_value}</td>"
            table_html += "</tr>"
        
        table_html += """
                    </tbody>
                </table>
            </div>
        </div>
        """
        
        return table_html
    
class PaginatedDataFrameTable(DataFrameTable):
    """A paginated version of DataFrameTable for handling large datasets."""
    
    def __init__(self, data, title, description, rows_per_page=20, **kwargs):
        super().__init__(data=data, title=title, description=description, **kwargs)
        self.rows_per_page = rows_per_page
        self.total_rows = len(data)
        self.total_pages = math.ceil(self.total_rows / self.rows_per_page)
        self.table_id = f"table-{hash(title)}-{id(self)}"
    
    def generate(self):
        from IPython.display import HTML
        
        if self.total_rows <= self.rows_per_page:
            # If data fits in one page, use the regular table
            return super().generate()
        
        # Generate paginated table
        table_html = f"""
        <div class="table-section">
            <h4>{self.title}</h4>
            <p class="text-muted">{self.description}</p>
            
            <!-- Pagination info -->
            <div class="d-flex justify-content-between align-items-center mb-2">
                <small class="text-muted">
                    Showing <span id="{self.table_id}-start">1</span> to 
                    <span id="{self.table_id}-end">{min(self.rows_per_page, self.total_rows)}</span> 
                    of {self.total_rows} entries
                </small>
                <div class="pagination-controls">
                    <button class="btn btn-sm btn-outline-secondary" id="{self.table_id}-prev" disabled>Previous</button>
                    <span class="mx-2">Page <span id="{self.table_id}-current">1</span> of {self.total_pages}</span>
                    <button class="btn btn-sm btn-outline-secondary" id="{self.table_id}-next" {'disabled' if self.total_pages <= 1 else ''}>Next</button>
                </div>
            </div>
            
            <div class="{self.style.get_container_class()}">
                <table class="{self.style.get_table_class()}" id="{self.table_id}">
                    <thead class="{self.style.get_header_class()}">
                        <tr>
        """
        
        # Add sortable headers
        for col_name, col in self.columns.items():
            table_html += f'''<th scope="col" title="{col.description}" class="sortable-header" data-column="{col_name}" style="cursor: pointer; user-select: none;">
                {col.name} <span class="sort-indicator"></span>
            </th>'''
        
        table_html += """
                        </tr>
                    </thead>
                    <tbody id="{}-body">
        """.format(self.table_id)
        
        # Add all data rows (we'll show/hide via JavaScript)
        display_data = self.data.reset_index(drop=True)
        for idx, (_, row) in enumerate(display_data.iterrows()):
            page_num = idx // self.rows_per_page + 1
            row_class = self.style.get_row_class()
            table_html += f'<tr class="{row_class}" data-page="{page_num}" style="{"" if page_num == 1 else "display: none;"}">'
            for col_name, col in self.columns.items():
                formatted_value = col.format_value(row[col_name])
                table_html += f"<td>{formatted_value}</td>"
            table_html += "</tr>"
        
        table_html += f"""
                    </tbody>
                </table>
            </div>
        </div>
        
        <style>
        .sortable-header:hover {{
            background-color: #f8f9fa;
        }}
        
        .sort-indicator {{
            font-size: 0.8em;
            color: #6c757d;
            margin-left: 5px;
        }}
        
        .sort-asc .sort-indicator::after {{
            content: '▲';
        }}
        
        .sort-desc .sort-indicator::after {{
            content: '▼';
        }}
        </style>
        
        <script>
        (function() {{
            const tableId = '{self.table_id}';
            const rowsPerPage = {self.rows_per_page};
            const totalPages = {self.total_pages};
            let currentPage = 1;
            let currentSort = {{ column: null, direction: 'asc' }};
            
            const prevBtn = document.getElementById(tableId + '-prev');
            const nextBtn = document.getElementById(tableId + '-next');
            const currentSpan = document.getElementById(tableId + '-current');
            const startSpan = document.getElementById(tableId + '-start');
            const endSpan = document.getElementById(tableId + '-end');
            const tbody = document.getElementById(tableId + '-body');
            const table = document.getElementById(tableId);
            
            // Store original row data for sorting
            const originalRows = Array.from(tbody.querySelectorAll('tr')).map((row, idx) => ({{
                element: row.cloneNode(true),
                data: {{}},
                originalIndex: idx
            }}));
            
            // Extract data from each row for sorting
            const headers = table.querySelectorAll('th[data-column]');
            headers.forEach((header, colIdx) => {{
                const columnName = header.dataset.column;
                originalRows.forEach((rowData, rowIdx) => {{
                    const cellText = rowData.element.cells[colIdx].textContent.trim();
                    // Try to parse as number, otherwise keep as string
                    rowData.data[columnName] = isNaN(cellText) || cellText === '' ? cellText : parseFloat(cellText);
                }});
            }});
            
            function sortData(column, direction) {{
                return [...originalRows].sort((a, b) => {{
                    let aVal = a.data[column];
                    let bVal = b.data[column];
                    
                    // Handle different data types
                    if (typeof aVal === 'string' && typeof bVal === 'string') {{
                        aVal = aVal.toLowerCase();
                        bVal = bVal.toLowerCase();
                    }}
                    
                    let comparison = 0;
                    if (aVal < bVal) comparison = -1;
                    else if (aVal > bVal) comparison = 1;
                    
                    return direction === 'desc' ? -comparison : comparison;
                }});
            }}
            
            function redistributePages(sortedRows) {{
                // Clear current tbody
                tbody.innerHTML = '';
                
                // Add sorted rows back with new page assignments
                sortedRows.forEach((rowData, idx) => {{
                    const newPageNum = Math.floor(idx / rowsPerPage) + 1;
                    const row = rowData.element.cloneNode(true);
                    row.setAttribute('data-page', newPageNum);
                    row.style.display = newPageNum === currentPage ? '' : 'none';
                    tbody.appendChild(row);
                }});
                
                // Update total pages if needed
                const newTotalPages = Math.ceil(sortedRows.length / rowsPerPage);
                if (currentPage > newTotalPages) {{
                    currentPage = Math.max(1, newTotalPages);
                }}
            }}
            
            function updateSortIndicators(activeColumn, direction) {{
                // Clear all sort indicators
                headers.forEach(header => {{
                    header.classList.remove('sort-asc', 'sort-desc');
                }});
                
                // Set active sort indicator
                if (activeColumn) {{
                    const activeHeader = table.querySelector(`th[data-column="${{activeColumn}}"]`);
                    if (activeHeader) {{
                        activeHeader.classList.add(direction === 'asc' ? 'sort-asc' : 'sort-desc');
                    }}
                }}
            }}
            
            function updateTable() {{
                // Hide all rows
                const rows = tbody.querySelectorAll('tr');
                rows.forEach(row => row.style.display = 'none');
                
                // Show current page rows
                const pageRows = tbody.querySelectorAll(`tr[data-page="${{currentPage}}"]`);
                pageRows.forEach(row => row.style.display = '');
                
                // Update pagination info
                const totalCurrentRows = tbody.querySelectorAll('tr').length;
                const start = (currentPage - 1) * rowsPerPage + 1;
                const end = Math.min(currentPage * rowsPerPage, totalCurrentRows);
                const currentTotalPages = Math.ceil(totalCurrentRows / rowsPerPage);
                
                startSpan.textContent = start;
                endSpan.textContent = end;
                currentSpan.textContent = currentPage;
                
                // Update button states
                prevBtn.disabled = currentPage === 1;
                nextBtn.disabled = currentPage === currentTotalPages || currentTotalPages === 0;
            }}
            
            // Add click event listeners to sortable headers
            headers.forEach(header => {{
                header.addEventListener('click', () => {{
                    const column = header.dataset.column;
                    
                    // Determine sort direction
                    if (currentSort.column === column) {{
                        currentSort.direction = currentSort.direction === 'asc' ? 'desc' : 'asc';
                    }} else {{
                        currentSort.column = column;
                        currentSort.direction = 'asc';
                    }}
                    
                    // Sort the data
                    const sortedRows = sortData(currentSort.column, currentSort.direction);
                    
                    // Reset to first page after sorting
                    currentPage = 1;
                    
                    // Redistribute rows across pages
                    redistributePages(sortedRows);
                    
                    // Update sort indicators
                    updateSortIndicators(currentSort.column, currentSort.direction);
                    
                    // Update table display
                    updateTable();
                }});
            }});
            
            prevBtn.addEventListener('click', () => {{
                if (currentPage > 1) {{
                    currentPage--;
                    updateTable();
                }}
            }});
            
            nextBtn.addEventListener('click', () => {{
                const totalCurrentRows = tbody.querySelectorAll('tr').length;
                const currentTotalPages = Math.ceil(totalCurrentRows / rowsPerPage);
                if (currentPage < currentTotalPages) {{
                    currentPage++;
                    updateTable();
                }}
            }});
            
            // Initial table update
            updateTable();
        }})();
        </script>
        """
        
        return table_html

class TabVisualizer:
    """Manages tabbed visualization of plot groups with simple, clean styling."""
    def __init__(self, title: str = "Spacemake QC Report"):
        self.title = title
        self.plot_groups: List[PlotGroup] = []
        self.unique_id = str(uuid.uuid4())[:8]
    
    def add_plot_group(self, group: PlotGroup) -> None:
        """Add a plot group to the visualizer."""
        self.plot_groups.append(group)
    
    def generate_html(self):
        """Generate complete HTML with simple, clean styling."""
        from IPython.display import HTML

        tabs_html = f'<ul class="nav nav-tabs nav-tabs-{self.unique_id}" role="tablist">'
        
        content_html = f'<div class="tab-content tab-content-{self.unique_id} mt-3">'
        
        for i, group in enumerate(self.plot_groups):
            active = i == 0
            clean_id = group.get_clean_id()
            unique_clean_id = f"{clean_id}-{self.unique_id}"
            
            # Add tab navigation item
            tabs_html += f"""
            <li class="nav-item" role="presentation">
                <button class="nav-link {'active' if active else ''}"
                        id="{unique_clean_id}-tab"
                        data-bs-toggle="tab"
                        data-bs-target="#{unique_clean_id}"
                        data-tab-container="{self.unique_id}"
                        type="button"
                        role="tab"
                        aria-controls="{unique_clean_id}"
                        aria-selected="{'true' if active else 'false'}">
                    {group.name}
                </button>
            </li>
            """
            
            content_html += f"""
            <div class="tab-pane fade {'show active' if active else ''}" 
                 id="{unique_clean_id}" 
                 role="tabpanel" 
                 aria-labelledby="{unique_clean_id}-tab">
                <h3>{group.name}</h3>
                <p class="lead text-muted">{group.description}</p>
                {chr(10).join(plot.generate() for plot in group.plots)}
            </div>
            """
        
        tabs_html += '</ul>'
        content_html += '</div>'

        html_template = f"""
        <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/css/bootstrap.min.css" rel="stylesheet">
        
        <div class="container-fluid tabvisualizer-{self.unique_id}">
            <div class="row">
                <div class="col-12">
                    <div class="card-body">
                        {tabs_html}
                        {content_html}
                    </div>
                </div>
            </div>
        </div>
        
        <style>
        .plot-section {{
            background: #f8f9fa;
            border-radius: 0.5rem;
            padding: 1.5rem;
            margin-bottom: 1.5rem;
            border: 1px solid #dee2e6;
        }}
        
        .table-section {{
            background: #f8f9fa;
            border-radius: 0.5rem;
            padding: 1.5rem;
            margin-bottom: 1.5rem;
            border: 1px solid #dee2e6;
        }}
        
        .table-container {{
            overflow-x: auto;
        }}
        
        .nav-tabs-{self.unique_id} .nav-link {{
            color: #495057;
            cursor: pointer;
        }}

        .nav-tabs-{self.unique_id} .nav-link::before {{
            content: none !important;
        }}

        .nav-tabs-{self.unique_id} .nav-item {{
            list-style: none !important;
        }}

        .nav-tabs-{self.unique_id} {{
            list-style-type: none !important;
        }}
        
        .nav-tabs-{self.unique_id} .nav-link.active {{
            color: #495057;
            background-color: #fff;
            border-color: #dee2e6 #dee2e6 #fff;
        }}
        
        .card {{
            border: none;
            box-shadow: 0 0.125rem 0.25rem rgba(0, 0, 0, 0.075);
        }}
        
        .table {{
            background: white;
        }}
        </style>
        
        <script src="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/js/bootstrap.bundle.min.js"></script>
        
        <script>
        // Initialize Bootstrap tabs manually for this specific TabVisualizer instance
        document.addEventListener('DOMContentLoaded', function() {{
            const tabButtons_{self.unique_id} = document.querySelectorAll('.tabvisualizer-{self.unique_id} [data-bs-toggle="tab"]');
            
            tabButtons_{self.unique_id}.forEach(button => {{
                button.addEventListener('click', function(e) {{
                    e.preventDefault();
                    
                    const containerId = this.getAttribute('data-tab-container');
                    const container = document.querySelector('.tabvisualizer-' + containerId);
                    
                    if (container) {{
                        // Remove active class from all tabs in THIS container only
                        container.querySelectorAll('.nav-tabs-{self.unique_id} .nav-link').forEach(link => {{
                            link.classList.remove('active');
                            link.setAttribute('aria-selected', 'false');
                        }});

                        // Hide tab content in THIS specific container only
                        container.querySelectorAll('.tab-content-{self.unique_id} .tab-pane').forEach(pane => {{
                            pane.classList.remove('show', 'active');
                        }});
                        
                        // Add active class to clicked tab
                        this.classList.add('active');
                        this.setAttribute('aria-selected', 'true');
                        
                        // Show corresponding content
                        const targetId = this.getAttribute('data-bs-target');
                        const targetPane = container.querySelector(targetId);
                        if (targetPane) {{
                            targetPane.classList.add('show', 'active');
                        }}
                    }}
                }});
            }});
        }});
        </script>
        """

        return HTML(html_template)

def histogram(values, axis, nbins=100, color="#000000", log=False, auto_log=True):
    # decide linear or logarithmic scale
    min_difference = values.max() - values.min()

    if log == False or (auto_log == True and np.abs(min_difference) < 100):
        hist, bins = np.histogram(values, bins=nbins)
        width = bins[1] - bins[0]
        axis.bar(bins[:-1], hist, width=width, color=color, align='edge')
    else:
        logbins = np.logspace(np.log10(bins[0] + 1), np.log10(bins[-1]), nbins)
        axis.hist(values, bins=logbins, color=color)
        axis.set_xscale("log")

    axis.spines[["right", "top"]].set_visible(False)

@ensure_anndata
def _scales_for_spatial_plot(adata):
    n_cells = len(adata)

    try:
        px_by_um = adata.uns["puck_variables"]["coord_by_um"]
        spot_diameter_um = adata.uns["puck_variables"]["spot_diameter_um"]
    except:
        logging.warning("Could not find 'coord_by_um' in the 'puck_variables' from the AnnData file. Setting to 1 as default")
        px_by_um = 1
        spot_diameter_um = 1
    
    meshed = False
    if "mesh_variables" in adata.uns.keys():
        meshed = True
        mesh_spot_diameter_um = adata.uns["mesh_variables"]["spot_diameter_um"]

    # Set limits and axes units for the spatial plots
    x_limits = adata.obsm["spatial"][:, 0].min(), adata.obsm["spatial"][:, 0].max()
    y_limits = adata.obsm["spatial"][:, 1].min(), adata.obsm["spatial"][:, 1].max()
    puck_width_um = (x_limits[1] - x_limits[0]) / px_by_um

    ratio = (x_limits[1] - x_limits[0]) / (y_limits[1] - y_limits[0])

    scale_factor = 2 if puck_width_um < 3000 else 3
    mm_dist = max(10**scale_factor, round(puck_width_um / 3, scale_factor))
    mm_diff = mm_dist / 1000

    def_plot_bead_size = 0.5 if n_cells > 5000 else 0.75
    def_plot_bead_size = 0.1 if n_cells > 10000 else def_plot_bead_size
    def_plot_bead_size = 0.05 if n_cells > 25000 else def_plot_bead_size

    puck_bead_size = max(
        def_plot_bead_size, mesh_spot_diameter_um if meshed else spot_diameter_um
    )
    x_mm_breaks = np.arange(0, puck_width_um, mm_dist)
    x_mm_breaks = [f"{round(i, 1)} mm" for i in x_mm_breaks * mm_diff / mm_dist]
    y_mm_breaks = np.arange(0, puck_width_um / ratio, mm_dist)
    y_mm_breaks = [f"{round(i, 1)} mm" for i in y_mm_breaks * mm_diff / mm_dist]

    x_breaks = np.arange(x_limits[0], x_limits[1], px_by_um * mm_dist)
    y_breaks = np.arange(y_limits[0], y_limits[1], px_by_um * mm_dist)

    return {"x_mm_breaks": x_mm_breaks,
            "y_mm_breaks": y_mm_breaks,
            "x_breaks": x_breaks,
            "y_breaks": y_breaks,
            "puck_bead_size": puck_bead_size,
            "px_by_um": px_by_um}

@ensure_anndata
def spatial(adata, spot_size=1.5, color="total_counts", cmap="magma", figsize=(5, 5), return_fig=True, **kwargs):
    if spot_size <= 0:
        raise ValueError("spot_size must be > 0")
    
    if (color not in adata.obs) and (color not in adata.var) and (color not in adata.var_names):
        logger.warning(f"No '{color}' found in adata")
        return None
    
    fig, axes = plt.subplots(1, 1, figsize=figsize)
    if len(adata) < 2:
        return fig, axes

    scales = _scales_for_spatial_plot(adata)

    sc.pl.spatial(
        adata,
        img_key=None,
        size=scales["puck_bead_size"] * spot_size,
        spot_size=scales["px_by_um"],
        color=color,
        ax=axes,
        show=False,
        cmap=cmap,
        **kwargs
    )
    axes.spines[["right", "top"]].set_visible(False)
    axes.set_xticks(scales["x_breaks"])
    axes.set_xticklabels(scales["x_mm_breaks"])
    axes.set_yticks(scales["y_breaks"])
    axes.set_yticklabels(scales["y_mm_breaks"])
    axes.set_ylabel("")
    axes.set_xlabel("")

    if return_fig:
        return fig, axes

@ensure_anndata
def neighborhood_enrichment(adata, key=None, spot_size=1.5, color="total_counts", cmap="magma", figsize=(4, 4), return_fig=True, **kwargs):
    if key is None or key == "":
        raise ValueError("`key` must be a valid .uns key")
    
    nhood_enrich_current_res = pd.DataFrame(adata.uns[f'{key}']['zscore'])
    nhood_enrich_current_res = pd.melt(nhood_enrich_current_res.reset_index(), id_vars='index')\
        .rename(columns={'index': 'cluster_a',
                         'variable': 'cluster_b',
                         'value': 'zscore'})

    ne_data = np.zeros(
        (
            nhood_enrich_current_res["cluster_a"].max() + 1,
            nhood_enrich_current_res["cluster_b"].max() + 1,
        )
    )
    ne_data[
        nhood_enrich_current_res["cluster_a"],
        nhood_enrich_current_res["cluster_b"],
    ] = nhood_enrich_current_res["zscore"]
    
    fig, axes = plt.subplots(1, 1, figsize=figsize)
    plmat = axes.matshow(
        ne_data, cmap="magma", vmin=-50, vmax=100, origin="lower"
    )
    cbar = plt.colorbar(plmat, fraction=0.046)
    cbar.set_label("Neighbor enrichment score")
    axes.set_xlabel("cluster identity")
    axes.set_ylabel("cluster identity")
    axes.spines[["right", "top"]].set_visible(False)
    axes.xaxis.tick_bottom()
    plt.tight_layout()

    if return_fig:
        return fig, axes

@ensure_anndata
def umap(adata, color=None, return_fig=True, **kwargs):
    if color is None or color == "":
        raise ValueError("`color` must be a valid adata.obs or adata.var key")

    fig, axes = plt.subplots(1, 1, figsize=(8, 5))
    sc.pl.umap(adata, color=f"{color}", ax=axes, show=False, **kwargs)
    axes.spines[["right", "top", "left", "bottom"]].set_visible(False)
    axes.set_ylabel("UMAP 0")
    axes.set_xlabel("UMAP 1")
    plt.tight_layout()

    if return_fig:
        return fig, axes

@ensure_anndata
def marker_gene_table(adata, rank_key=None):
    if rank_key is None or rank_key == "":
        raise ValueError("`rank_key` must be a valid adata.uns key")

    if not 'names' in adata.uns[rank_key]:
        return None

    df = pd.DataFrame(adata.uns[rank_key]['names'])\
            .melt(var_name = 'cluster', value_name = 'gene')

    for key in ['logfoldchanges', 'pvals', 'pvals_adj']:
        df_key = pd.DataFrame(adata.uns[rank_key][key])\
            .melt(var_name = 'cluster', value_name = key)
        df[key] = df_key[key]
        
    df.set_index(['gene', 'cluster'], inplace=True)
        
    for key in ['pts', 'pts_rest']:
        df2 = adata.uns[rank_key][key]
        df2['gene'] = df2.index
        df2 = df2.melt(var_name='cluster', id_vars='gene')\
            .set_index(['gene', 'cluster'])
        
        df[key] = df2.loc[df.index].value
            
    df.reset_index(inplace=True)
    return df

@ensure_anndata
def knee_plot(adata, figsize=(5, 3), return_fig=True):
    fig, axis = plt.subplots(1, 1, figsize=figsize)
    axis.plot(
        np.arange(len(adata)),
        np.cumsum(adata.obs["n_reads"].values),
        color="black",
        linewidth=1,
    )
    axis.set_ylabel("Cumulative\nsum of reads")
    axis.set_xlabel("Beads sorted by number of reads")
    axis.spines[["right", "top"]].set_visible(False)
    plt.tight_layout()

    if return_fig:
        return fig, axis

@ensure_anndata
def umi_cutoff(adata, figsize=(8, 4), return_fig=True):
    umi_cutoffs = np.arange(10, 20000, 10)

    def summarise_dge_summary(df, umi_cutoff):
        df_filter = df[df["total_counts"] > umi_cutoff]

        df_summary = df_filter[
            [
                "n_reads",
                "total_counts",
                "n_genes_by_counts",
                "reads_per_counts",
            ]
        ].median()
        df_summary["n_beads"] = len(df_filter)

        return df_summary

    umi_cutoff_data = pd.DataFrame(
        np.vstack([summarise_dge_summary(adata.obs, umi_cutoff).values for umi_cutoff in umi_cutoffs]),
        columns=[
            "n_reads",
            "total_counts",
            "n_genes_by_counts",
            "reads_per_counts",
            "n_beads",
        ],
    )
    umi_cutoff_data["umi_cutoffs"] = umi_cutoffs

    fig, axes = plt.subplots(2, 3, figsize=figsize)
    axes[0, 0].plot(
        umi_cutoff_data["umi_cutoffs"],
        umi_cutoff_data["n_beads"],
        color="black",
        linewidth=1,
    )
    axes[0, 0].set_ylabel("number of\nspatial units")
    axes[0, 0].set_xlabel("minimum UMI")
    axes[0, 0].set_yscale("log")
    axes[0, 0].set_xscale("log")
    axes[0, 0].spines[["right", "top"]].set_visible(False)

    axes[0, 1].plot(
        umi_cutoff_data["umi_cutoffs"],
        umi_cutoff_data["n_reads"],
        color=metrics_colors["reads"],
        linewidth=1,
    )
    axes[0, 1].set_ylabel("median reads\nper spatial unit")
    axes[0, 1].set_xlabel("minimum UMI")
    axes[0, 1].set_xscale("log")
    axes[0, 1].spines[["right", "top"]].set_visible(False)

    axes[0, 2].plot(
        umi_cutoff_data["umi_cutoffs"],
        umi_cutoff_data["n_genes_by_counts"],
        color=metrics_colors["genes"],
        linewidth=1,
    )
    axes[0, 2].set_ylabel("median genes\nper spatial unit")
    axes[0, 2].set_xlabel("minimum UMI")
    axes[0, 2].set_xscale("log")
    axes[0, 2].spines[["right", "top"]].set_visible(False)

    axes[1, 0].plot(
        umi_cutoff_data["umi_cutoffs"],
        umi_cutoff_data["total_counts"],
        color=metrics_colors["umis"],
        linewidth=1,
    )
    axes[1, 0].set_ylabel("median UMIs\nper spatial unit")
    axes[1, 0].set_xlabel("minimum UMI")
    axes[1, 0].set_xscale("log")
    axes[1, 0].spines[["right", "top"]].set_visible(False)

    axes[1, 1].plot(
        umi_cutoff_data["umi_cutoffs"],
        umi_cutoff_data["reads_per_counts"],
        color=metrics_colors["pcr"],
        linewidth=1,
    )
    axes[1, 1].set_ylabel("median reads/UMI\nper spatial unit")
    axes[1, 1].set_xlabel("minimum UMI")
    axes[1, 1].set_xscale("log")
    axes[1, 1].spines[["right", "top"]].set_visible(False)
    axes[1, 2].axis("off")

    plt.tight_layout()

    if return_fig:
        return fig, axes

@ensure_anndata
def plot_histogram_beads(adata, figsize=(7, 3.5), return_fig=True):
    fig, axes = plt.subplots(2, 2, figsize=figsize)
    histogram(adata.obs["n_reads"], axes[0, 0], 100, metrics_colors["reads"])
    axes[0, 0].set_ylabel("# of\nspatial units")
    axes[0, 0].set_xlabel("# of reads")

    reads_per_counts = adata.obs["reads_per_counts"]
    reads_per_counts = np.nan_to_num(reads_per_counts)
    histogram(reads_per_counts, axes[0, 1], 100, metrics_colors["pcr"])
    axes[0, 1].set_ylabel("# of\nspatial units")
    axes[0, 1].set_xlabel("# of reads/UMI")

    histogram(adata.obs["n_genes_by_counts"], axes[1, 0], 100, metrics_colors["genes"])
    axes[1, 0].set_ylabel("# of\nspatial units")
    axes[1, 0].set_xlabel("# of genes")

    histogram(adata.obs["total_counts"], axes[1, 1], 100, metrics_colors["umis"])
    axes[1, 1].set_ylabel("# of\nspatial units")
    axes[1, 1].set_xlabel("# of UMIs")

    plt.tight_layout()

    if return_fig:
        return fig, axes

@ensure_anndata
def nucleotide_distribution_per_bead(adata, figsize=(8, 4), return_fig=True):
    adata.obs["reads_cumsum"] = adata.obs["n_reads"].cumsum()
    adata.obs["quartile"] = pd.cut(
        adata.obs["reads_cumsum"],
        bins=4,
        include_lowest=True,
        labels=["Q1", "Q2", "Q3", "Q4"],
    )
    adata.obs.drop(columns=["reads_cumsum"], inplace=True)

    # Create a dataframe to show the count of nucleotides/barcode
    cell_bc_len = len(adata.obs["cell_bc"].iloc[0])
    nucls = adata.obs["cell_bc"].str.strip().apply(list).apply(pd.Series)
    nucls = pd.concat([adata.obs[["cell_bc", "quartile"]], nucls], axis=1)
    nucls = nucls.melt(id_vars=["cell_bc", "quartile"], var_name="pos", value_name="nucl")
    nucls = nucls.groupby(["pos", "nucl", "quartile"]).size().reset_index(name="nucl_count")
    nucls = nucls.pivot_table(
        index=["pos", "nucl"], columns="quartile", values="nucl_count", fill_value=0
    ).reset_index()
    lbl_df = adata.obs.groupby("quartile").size().reset_index(name="lbl")
    lbl_df["lbl"] = lbl_df.apply(lambda row: f"{row['quartile']} (n={row['lbl']})", axis=1)
    lbls = dict(zip(lbl_df["quartile"], lbl_df["lbl"]))

    # Create the plot
    fig, axes = plt.subplots(2, 2, figsize=figsize)
    x = np.arange(1, cell_bc_len + 1)  # the label locations

    for name, group in nucls.groupby("nucl"):
        axes[0, 0].plot(x, group["Q1"], color=nucl_clrs[name], linewidth=2, label=name)
        axes[0, 0].set_xlim(0.1, cell_bc_len + 1.5)
        axes[0, 0].set_xticks([])
        axes[0, 0].set_title(lbls["Q1"])
        axes[0, 0].spines[["right", "top", "bottom"]].set_visible(False)

        axes[0, 1].plot(x, group["Q2"], color=nucl_clrs[name], linewidth=2)
        axes[0, 1].set_xlim(0.1, cell_bc_len + 1.5)
        axes[0, 1].set_xticks([])
        axes[0, 1].set_title(lbls["Q2"])
        axes[0, 1].spines[["right", "top", "bottom"]].set_visible(False)

        axes[1, 0].plot(x, group["Q3"], color=nucl_clrs[name], linewidth=2)
        axes[1, 0].set_xlim(0.1, cell_bc_len + 1.5)
        axes[1, 0].set_title(lbls["Q3"])
        axes[1, 0].spines[["right", "top"]].set_visible(False)
        axes[1, 0].set_xticks(list(set(list(range(1, cell_bc_len + 1, 2)) + [cell_bc_len])))

        axes[1, 1].plot(x, group["Q4"], color=nucl_clrs[name], linewidth=2)
        axes[1, 1].set_xlim(0.1, cell_bc_len + 1.5)
        axes[1, 1].set_title(lbls["Q4"])
        axes[1, 1].spines[["right", "top"]].set_visible(False)
        axes[1, 1].set_xticks(list(set(list(range(1, cell_bc_len + 1, 2)) + [cell_bc_len])))

    handles, labels = axes[0, 0].get_legend_handles_labels()
    legend = fig.legend(
        handles,
        labels,
        loc="center right",
        bbox_to_anchor=(1.1, 0.5),
        borderaxespad=0,
        title="Nucleotide",
    )
    legend.set_frame_on(False)
    fig.text(0.5, 0.0, "nucleotide position in the barcode", ha="center")
    plt.tight_layout()

    if return_fig:
        return fig, axes

@ensure_anndata
def entropy_compression(adata, nbins=30, figsize=(7, 4), return_fig=True):
    fig, axes = plt.subplots(2, 1, figsize=figsize)
    axes[0].hist(
        adata.obs["theoretical_entropy"].values,
        color=cpalette["grey"],
        bins=nbins,
        edgecolor="black",
        label="Theoretical",
    )
    axes[0].hist(
        adata.obs["exact_entropy"].values,
        color=cpalette["orange"],
        bins=nbins,
        edgecolor="black",
        alpha=0.7,
        label="Observed",
    )
    axes[0].set_xlabel("Shannon entropy of barcodes")
    axes[0].set_ylabel("# of barcodes")
    axes[0].spines[["right", "top"]].set_visible(False)
    legend0 = axes[0].legend(loc="upper left")
    legend0.set_frame_on(False)

    axes[1].hist(
        adata.obs["theoretical_compression"].values,
        color=cpalette["grey"],
        bins=nbins,
        edgecolor="black",
        label="Theoretical",
    )
    axes[1].hist(
        adata.obs["exact_compression"].values,
        color=cpalette["orange"],
        bins=nbins,
        edgecolor="black",
        alpha=0.7,
        label="Observed",
    )
    axes[1].set_xlabel("Length of barcodes after compression")
    axes[1].set_ylabel("# of barcodes")
    axes[1].spines[["right", "top"]].set_visible(False)
    axes[1].set_xlim(left=0)

    if return_fig:
        return fig, axes


def density_per_downsampling(values, metric, log_scale=True, color="#000000", title="", figsize=None, return_fig=True):
    if figsize is None:
        figsize = (5, 0.5 * len(PCT_DOWNSAMPLE_TO_PLOT))
        
    fig, axes = plt.subplots(len(PCT_DOWNSAMPLE_TO_PLOT), 1, figsize=figsize)

    i = 0
    for downsample_pct, value_density in values.groupby("_downsample_pct_report"):
        if int(downsample_pct) in PCT_DOWNSAMPLE_TO_PLOT:
            density_function = gaussian_kde(np.nan_to_num(value_density[metric]), bw_method=0.1)
            x = np.linspace(1, max(np.nan_to_num(values[metric])), 100)

            axes[i].plot(x, density_function(x), color="black", linewidth=1)
            axes[i].fill_between(x, density_function(x), color=color)
            axes[i].set_yticks([])

            if log_scale:
                axes[i].set_xscale("log")

            axes[i].spines[["right", "top", "bottom"]].set_visible(False)
            axes[i].text(1.05, 0.5, f"{downsample_pct}%", transform=axes[i].transAxes, va="center")
            i += 1

        axes[-1].spines[["right", "top"]].set_visible(False)
        axes[-1].spines[["left", "bottom"]].set_visible(True)
        axes[-1].set_xlabel(title)

    for i in range(i - 1):
        axes[i].set_xticks([])

    fig.text(0.0, 0.6, "density", va="center", rotation="vertical")
    plt.tight_layout()

    if return_fig:
        return fig, axes


def median_per_downsampling(values, metric, umi_cutoffs, color="#000000", title="", figsize=(5, 3), return_fig=True):
    fig, axes = plt.subplots(1, 1, figsize=figsize)

    lines = ["-", "--", "-.", ":"]
    linecycler = cycle(lines)
    handles, labels = [], []

    for umi_cutoff in umi_cutoffs:
        _values = values[values["total_counts"] > umi_cutoff]
        median_values = (
            _values[[metric, "_downsample_pct_report"]].groupby("_downsample_pct_report").median().reset_index()
        )

        linestyle = next(linecycler)

        (line,) = axes.plot(
            median_values["_downsample_pct_report"], median_values[metric], linestyle, color=color, label=umi_cutoff
        )
        axes.scatter(
            median_values["_downsample_pct_report"], median_values[metric], s=20, color=color, edgecolors="black"
        )

        handles.append(line)
        labels.append(umi_cutoff)

    axes.set_xticks(PCT_DOWNSAMPLE_TO_PLOT)
    axes.set_xticklabels([f"{pct}%" for pct in PCT_DOWNSAMPLE_TO_PLOT])
    axes.spines[["right", "top"]].set_visible(False)
    axes.set_xlabel("downsampling percentage")
    axes.set_ylabel(title)

    legend = axes.legend(handles, labels, loc="lower right", title="UMI cutoff")
    legend.set_frame_on(False)

    plt.tight_layout()
    if return_fig:
        return fig, axes


def deciled_median(decile_dat, figsize=(6, 4), return_fig=True):
    fig, axes = plt.subplots(3, 2, figsize=figsize)

    # Iterate through each unique 'observation' for facetting
    for i, (obs, data) in enumerate(decile_dat.groupby("observation")):
        for _obs, _data in data.groupby("decile"):
            axes.flatten()[i].plot(
                _data["_downsample_pct_report"], _data["value"], label=_obs, linewidth=0.6, color=parula_dict[_obs]
            )
            axes.flatten()[i].scatter(
                _data["_downsample_pct_report"], _data["value"], s=20, edgecolors="black", color=parula_dict[_obs]
            )

        axes.flatten()[i].set_xticks([0, 20, 40, 60, 80, 100])
        axes.flatten()[i].set_xticklabels(["0", "20", "40", "60", "80", "100"])
        axes.flatten()[i].set_title(obs)
        axes.flatten()[i].spines[["top", "right"]].set_visible(False)

    axes.flatten()[i].set_xlabel("downsampling percentage")
    if i % 2 == 0:
        axes.flatten()[-1].axis("off")

    # Create a single legend at the bottom
    handles, labels = [], []
    for obs in parula_dict:
        handles.append(plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=parula_dict[obs], markersize=8))
        labels.append(str(obs))

    fig.legend(handles, labels, title="Decile", loc="lower right", ncol=3, bbox_to_anchor=(0.95, 0.02))

    plt.tight_layout()

    if return_fig:
        return fig, axes