I/O and VOSA helper functions
=============================

Reading photometric data
------------------------

From VOSA files
~~~~~~~~~~~~~~~

.. code:: ipython3

    import src.binary_sed_fitting as bsf
    file_name = 'data/vosa_results/WOCS2002.bfit.phot.dat'
    data = bsf.load_data(file_name, mode='vosa')
    data


.. parsed-literal::

    /Users/vikrantjadhav/Documents/work/Binary_SED_Fitting/src/binary_sed_fitting.py:158: FutureWarning: The 'delim_whitespace' keyword in pd.read_csv is deprecated and will be removed in a future version. Use ``sep='\s+'`` instead
      data = pd.read_csv(file_name, engine='python', comment='#',
    



.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>wavelength</th>
          <th>flux</th>
          <th>error</th>
        </tr>
        <tr>
          <th>FilterID</th>
          <th></th>
          <th></th>
          <th></th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>Astrosat/UVIT.F148W</th>
          <td>1481.000000</td>
          <td>2.548070e-15</td>
          <td>1.361179e-16</td>
        </tr>
        <tr>
          <th>Astrosat/UVIT.F154W</th>
          <td>1541.000000</td>
          <td>2.285731e-15</td>
          <td>1.368403e-16</td>
        </tr>
        <tr>
          <th>Astrosat/UVIT.F169M</th>
          <td>1608.000000</td>
          <td>2.360782e-15</td>
          <td>1.261129e-16</td>
        </tr>
        <tr>
          <th>GALEX/GALEX.NUV</th>
          <td>2303.366368</td>
          <td>3.296702e-15</td>
          <td>8.677147e-17</td>
        </tr>
        <tr>
          <th>KPNO/Mosaic.B</th>
          <td>4357.276538</td>
          <td>9.309407e-14</td>
          <td>6.130611e-15</td>
        </tr>
        <tr>
          <th>GAIA/GAIA3.Gbp</th>
          <td>5035.750275</td>
          <td>1.028575e-13</td>
          <td>3.881303e-16</td>
        </tr>
        <tr>
          <th>KPNO/Mosaic.V</th>
          <td>5366.240786</td>
          <td>1.079375e-13</td>
          <td>7.108112e-15</td>
        </tr>
        <tr>
          <th>GAIA/GAIA3.G</th>
          <td>5822.388714</td>
          <td>9.174216e-14</td>
          <td>2.539155e-16</td>
        </tr>
        <tr>
          <th>GAIA/GAIA3.Grp</th>
          <td>7619.959993</td>
          <td>7.878734e-14</td>
          <td>3.154436e-16</td>
        </tr>
        <tr>
          <th>KPNO/Mosaic.I</th>
          <td>8101.609574</td>
          <td>4.470797e-14</td>
          <td>2.944196e-15</td>
        </tr>
        <tr>
          <th>GAIA/GAIA3.Grvs</th>
          <td>8578.159519</td>
          <td>6.896280e-14</td>
          <td>5.559650e-16</td>
        </tr>
        <tr>
          <th>2MASS/2MASS.J</th>
          <td>12350.000000</td>
          <td>3.612486e-14</td>
          <td>7.319890e-16</td>
        </tr>
        <tr>
          <th>2MASS/2MASS.H</th>
          <td>16620.000000</td>
          <td>1.882778e-14</td>
          <td>3.468206e-16</td>
        </tr>
        <tr>
          <th>2MASS/2MASS.Ks</th>
          <td>21590.000000</td>
          <td>7.827519e-15</td>
          <td>1.297694e-16</td>
        </tr>
        <tr>
          <th>WISE/WISE.W1</th>
          <td>33526.000000</td>
          <td>1.546087e-15</td>
          <td>2.990397e-17</td>
        </tr>
        <tr>
          <th>WISE/WISE.W2</th>
          <td>46028.000000</td>
          <td>4.313308e-16</td>
          <td>7.150866e-18</td>
        </tr>
        <tr>
          <th>WISE/WISE.W3</th>
          <td>115608.000000</td>
          <td>1.230361e-17</td>
          <td>4.306178e-19</td>
        </tr>
      </tbody>
    </table>
    </div>



From a CSV file
~~~~~~~~~~~~~~~

-  The flux has to be extinction corrected
-  units of flux = erg/s/A/cm2
-  units of wavelength = A

.. code:: ipython3

    import src.binary_sed_fitting as bsf
    file_name = 'data/extinction_corrected_flux_files/WOCS2002.csv'
    data = bsf.load_data(file_name, mode='csv')
    data




.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>wavelength</th>
          <th>flux</th>
          <th>error</th>
        </tr>
        <tr>
          <th>FilterID</th>
          <th></th>
          <th></th>
          <th></th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>Astrosat/UVIT.F148W</th>
          <td>1481.000000</td>
          <td>2.548070e-15</td>
          <td>1.361179e-16</td>
        </tr>
        <tr>
          <th>Astrosat/UVIT.F154W</th>
          <td>1541.000000</td>
          <td>2.285731e-15</td>
          <td>1.368403e-16</td>
        </tr>
        <tr>
          <th>Astrosat/UVIT.F169M</th>
          <td>1608.000000</td>
          <td>2.360782e-15</td>
          <td>1.261129e-16</td>
        </tr>
        <tr>
          <th>GALEX/GALEX.NUV</th>
          <td>2303.366368</td>
          <td>3.296702e-15</td>
          <td>8.677147e-17</td>
        </tr>
        <tr>
          <th>KPNO/Mosaic.B</th>
          <td>4357.276538</td>
          <td>9.309407e-14</td>
          <td>6.130611e-15</td>
        </tr>
        <tr>
          <th>GAIA/GAIA3.Gbp</th>
          <td>5035.750275</td>
          <td>1.028575e-13</td>
          <td>3.881303e-16</td>
        </tr>
        <tr>
          <th>KPNO/Mosaic.V</th>
          <td>5366.240786</td>
          <td>1.079375e-13</td>
          <td>7.108112e-15</td>
        </tr>
        <tr>
          <th>GAIA/GAIA3.G</th>
          <td>5822.388714</td>
          <td>9.174216e-14</td>
          <td>2.539155e-16</td>
        </tr>
        <tr>
          <th>GAIA/GAIA3.Grp</th>
          <td>7619.959993</td>
          <td>7.878734e-14</td>
          <td>3.154436e-16</td>
        </tr>
        <tr>
          <th>KPNO/Mosaic.I</th>
          <td>8101.609574</td>
          <td>4.470797e-14</td>
          <td>2.944196e-15</td>
        </tr>
        <tr>
          <th>GAIA/GAIA3.Grvs</th>
          <td>8578.159519</td>
          <td>6.896280e-14</td>
          <td>5.559650e-16</td>
        </tr>
        <tr>
          <th>2MASS/2MASS.J</th>
          <td>12350.000000</td>
          <td>3.612486e-14</td>
          <td>7.319890e-16</td>
        </tr>
        <tr>
          <th>2MASS/2MASS.H</th>
          <td>16620.000000</td>
          <td>1.882778e-14</td>
          <td>3.468206e-16</td>
        </tr>
        <tr>
          <th>2MASS/2MASS.Ks</th>
          <td>21590.000000</td>
          <td>7.827519e-15</td>
          <td>1.297694e-16</td>
        </tr>
        <tr>
          <th>WISE/WISE.W1</th>
          <td>33526.000000</td>
          <td>1.546087e-15</td>
          <td>2.990397e-17</td>
        </tr>
        <tr>
          <th>WISE/WISE.W2</th>
          <td>46028.000000</td>
          <td>4.313308e-16</td>
          <td>7.150866e-18</td>
        </tr>
        <tr>
          <th>WISE/WISE.W3</th>
          <td>115608.000000</td>
          <td>1.230361e-17</td>
          <td>4.306178e-19</td>
        </tr>
      </tbody>
    </table>
    </div>



Reading model files
-------------------

.. code:: ipython3

    import xarray as xr
    import src.binary_sed_fitting as bsf
    
    da = xr.open_dataarray(bsf.DIR_MODELS+'koester_synthetic_photometry.nc')
    da




.. raw:: html

    <div><svg style="position: absolute; width: 0; height: 0; overflow: hidden">
    <defs>
    <symbol id="icon-database" viewBox="0 0 32 32">
    <path d="M16 0c-8.837 0-16 2.239-16 5v4c0 2.761 7.163 5 16 5s16-2.239 16-5v-4c0-2.761-7.163-5-16-5z"></path>
    <path d="M16 17c-8.837 0-16-2.239-16-5v6c0 2.761 7.163 5 16 5s16-2.239 16-5v-6c0 2.761-7.163 5-16 5z"></path>
    <path d="M16 26c-8.837 0-16-2.239-16-5v6c0 2.761 7.163 5 16 5s16-2.239 16-5v-6c0 2.761-7.163 5-16 5z"></path>
    </symbol>
    <symbol id="icon-file-text2" viewBox="0 0 32 32">
    <path d="M28.681 7.159c-0.694-0.947-1.662-2.053-2.724-3.116s-2.169-2.030-3.116-2.724c-1.612-1.182-2.393-1.319-2.841-1.319h-15.5c-1.378 0-2.5 1.121-2.5 2.5v27c0 1.378 1.122 2.5 2.5 2.5h23c1.378 0 2.5-1.122 2.5-2.5v-19.5c0-0.448-0.137-1.23-1.319-2.841zM24.543 5.457c0.959 0.959 1.712 1.825 2.268 2.543h-4.811v-4.811c0.718 0.556 1.584 1.309 2.543 2.268zM28 29.5c0 0.271-0.229 0.5-0.5 0.5h-23c-0.271 0-0.5-0.229-0.5-0.5v-27c0-0.271 0.229-0.5 0.5-0.5 0 0 15.499-0 15.5 0v7c0 0.552 0.448 1 1 1h7v19.5z"></path>
    <path d="M23 26h-14c-0.552 0-1-0.448-1-1s0.448-1 1-1h14c0.552 0 1 0.448 1 1s-0.448 1-1 1z"></path>
    <path d="M23 22h-14c-0.552 0-1-0.448-1-1s0.448-1 1-1h14c0.552 0 1 0.448 1 1s-0.448 1-1 1z"></path>
    <path d="M23 18h-14c-0.552 0-1-0.448-1-1s0.448-1 1-1h14c0.552 0 1 0.448 1 1s-0.448 1-1 1z"></path>
    </symbol>
    </defs>
    </svg>
    <style>/* CSS stylesheet for displaying xarray objects in jupyterlab.
     *
     */
    
    :root {
      --xr-font-color0: var(--jp-content-font-color0, rgba(0, 0, 0, 1));
      --xr-font-color2: var(--jp-content-font-color2, rgba(0, 0, 0, 0.54));
      --xr-font-color3: var(--jp-content-font-color3, rgba(0, 0, 0, 0.38));
      --xr-border-color: var(--jp-border-color2, #e0e0e0);
      --xr-disabled-color: var(--jp-layout-color3, #bdbdbd);
      --xr-background-color: var(--jp-layout-color0, white);
      --xr-background-color-row-even: var(--jp-layout-color1, white);
      --xr-background-color-row-odd: var(--jp-layout-color2, #eeeeee);
    }
    
    html[theme=dark],
    body[data-theme=dark],
    body.vscode-dark {
      --xr-font-color0: rgba(255, 255, 255, 1);
      --xr-font-color2: rgba(255, 255, 255, 0.54);
      --xr-font-color3: rgba(255, 255, 255, 0.38);
      --xr-border-color: #1F1F1F;
      --xr-disabled-color: #515151;
      --xr-background-color: #111111;
      --xr-background-color-row-even: #111111;
      --xr-background-color-row-odd: #313131;
    }
    
    .xr-wrap {
      display: block !important;
      min-width: 300px;
      max-width: 700px;
    }
    
    .xr-text-repr-fallback {
      /* fallback to plain text repr when CSS is not injected (untrusted notebook) */
      display: none;
    }
    
    .xr-header {
      padding-top: 6px;
      padding-bottom: 6px;
      margin-bottom: 4px;
      border-bottom: solid 1px var(--xr-border-color);
    }
    
    .xr-header > div,
    .xr-header > ul {
      display: inline;
      margin-top: 0;
      margin-bottom: 0;
    }
    
    .xr-obj-type,
    .xr-array-name {
      margin-left: 2px;
      margin-right: 10px;
    }
    
    .xr-obj-type {
      color: var(--xr-font-color2);
    }
    
    .xr-sections {
      padding-left: 0 !important;
      display: grid;
      grid-template-columns: 150px auto auto 1fr 20px 20px;
    }
    
    .xr-section-item {
      display: contents;
    }
    
    .xr-section-item input {
      display: none;
    }
    
    .xr-section-item input + label {
      color: var(--xr-disabled-color);
    }
    
    .xr-section-item input:enabled + label {
      cursor: pointer;
      color: var(--xr-font-color2);
    }
    
    .xr-section-item input:enabled + label:hover {
      color: var(--xr-font-color0);
    }
    
    .xr-section-summary {
      grid-column: 1;
      color: var(--xr-font-color2);
      font-weight: 500;
    }
    
    .xr-section-summary > span {
      display: inline-block;
      padding-left: 0.5em;
    }
    
    .xr-section-summary-in:disabled + label {
      color: var(--xr-font-color2);
    }
    
    .xr-section-summary-in + label:before {
      display: inline-block;
      content: '►';
      font-size: 11px;
      width: 15px;
      text-align: center;
    }
    
    .xr-section-summary-in:disabled + label:before {
      color: var(--xr-disabled-color);
    }
    
    .xr-section-summary-in:checked + label:before {
      content: '▼';
    }
    
    .xr-section-summary-in:checked + label > span {
      display: none;
    }
    
    .xr-section-summary,
    .xr-section-inline-details {
      padding-top: 4px;
      padding-bottom: 4px;
    }
    
    .xr-section-inline-details {
      grid-column: 2 / -1;
    }
    
    .xr-section-details {
      display: none;
      grid-column: 1 / -1;
      margin-bottom: 5px;
    }
    
    .xr-section-summary-in:checked ~ .xr-section-details {
      display: contents;
    }
    
    .xr-array-wrap {
      grid-column: 1 / -1;
      display: grid;
      grid-template-columns: 20px auto;
    }
    
    .xr-array-wrap > label {
      grid-column: 1;
      vertical-align: top;
    }
    
    .xr-preview {
      color: var(--xr-font-color3);
    }
    
    .xr-array-preview,
    .xr-array-data {
      padding: 0 5px !important;
      grid-column: 2;
    }
    
    .xr-array-data,
    .xr-array-in:checked ~ .xr-array-preview {
      display: none;
    }
    
    .xr-array-in:checked ~ .xr-array-data,
    .xr-array-preview {
      display: inline-block;
    }
    
    .xr-dim-list {
      display: inline-block !important;
      list-style: none;
      padding: 0 !important;
      margin: 0;
    }
    
    .xr-dim-list li {
      display: inline-block;
      padding: 0;
      margin: 0;
    }
    
    .xr-dim-list:before {
      content: '(';
    }
    
    .xr-dim-list:after {
      content: ')';
    }
    
    .xr-dim-list li:not(:last-child):after {
      content: ',';
      padding-right: 5px;
    }
    
    .xr-has-index {
      font-weight: bold;
    }
    
    .xr-var-list,
    .xr-var-item {
      display: contents;
    }
    
    .xr-var-item > div,
    .xr-var-item label,
    .xr-var-item > .xr-var-name span {
      background-color: var(--xr-background-color-row-even);
      margin-bottom: 0;
    }
    
    .xr-var-item > .xr-var-name:hover span {
      padding-right: 5px;
    }
    
    .xr-var-list > li:nth-child(odd) > div,
    .xr-var-list > li:nth-child(odd) > label,
    .xr-var-list > li:nth-child(odd) > .xr-var-name span {
      background-color: var(--xr-background-color-row-odd);
    }
    
    .xr-var-name {
      grid-column: 1;
    }
    
    .xr-var-dims {
      grid-column: 2;
    }
    
    .xr-var-dtype {
      grid-column: 3;
      text-align: right;
      color: var(--xr-font-color2);
    }
    
    .xr-var-preview {
      grid-column: 4;
    }
    
    .xr-index-preview {
      grid-column: 2 / 5;
      color: var(--xr-font-color2);
    }
    
    .xr-var-name,
    .xr-var-dims,
    .xr-var-dtype,
    .xr-preview,
    .xr-attrs dt {
      white-space: nowrap;
      overflow: hidden;
      text-overflow: ellipsis;
      padding-right: 10px;
    }
    
    .xr-var-name:hover,
    .xr-var-dims:hover,
    .xr-var-dtype:hover,
    .xr-attrs dt:hover {
      overflow: visible;
      width: auto;
      z-index: 1;
    }
    
    .xr-var-attrs,
    .xr-var-data,
    .xr-index-data {
      display: none;
      background-color: var(--xr-background-color) !important;
      padding-bottom: 5px !important;
    }
    
    .xr-var-attrs-in:checked ~ .xr-var-attrs,
    .xr-var-data-in:checked ~ .xr-var-data,
    .xr-index-data-in:checked ~ .xr-index-data {
      display: block;
    }
    
    .xr-var-data > table {
      float: right;
    }
    
    .xr-var-name span,
    .xr-var-data,
    .xr-index-name div,
    .xr-index-data,
    .xr-attrs {
      padding-left: 25px !important;
    }
    
    .xr-attrs,
    .xr-var-attrs,
    .xr-var-data,
    .xr-index-data {
      grid-column: 1 / -1;
    }
    
    dl.xr-attrs {
      padding: 0;
      margin: 0;
      display: grid;
      grid-template-columns: 125px auto;
    }
    
    .xr-attrs dt,
    .xr-attrs dd {
      padding: 0;
      margin: 0;
      float: left;
      padding-right: 10px;
      width: auto;
    }
    
    .xr-attrs dt {
      font-weight: normal;
      grid-column: 1;
    }
    
    .xr-attrs dt:hover span {
      display: inline-block;
      background: var(--xr-background-color);
      padding-right: 10px;
    }
    
    .xr-attrs dd {
      grid-column: 2;
      white-space: pre-wrap;
      word-break: break-all;
    }
    
    .xr-icon-database,
    .xr-icon-file-text2,
    .xr-no-icon {
      display: inline-block;
      vertical-align: middle;
      width: 1em;
      height: 1.5em !important;
      stroke-width: 0;
      stroke: currentColor;
      fill: currentColor;
    }
    </style><pre class='xr-text-repr-fallback'>&lt;xarray.DataArray (FilterID: 8162, Te: 82, logg: 13)&gt; Size: 70MB
    [8700692 values with dtype=float64]
    Coordinates:
      * FilterID  (FilterID) &lt;U38 1MB &#x27;IUE/IUE.1250-1300&#x27; ... &#x27;QUIJOTE/MFI.11GHz_H3&#x27;
      * Te        (Te) int32 328B 5000 5250 5500 5750 ... 50000 60000 70000 80000
      * logg      (logg) float64 104B 6.5 6.75 7.0 7.25 7.5 ... 8.75 9.0 9.25 9.5
    Attributes:
        Wavelengths:  [1.28470318e+03 1.35336435e+03 1.35623969e+03 ... 3.8683349...
        unit:         erg/s/cm2/A
        long_name:    Flux</pre><div class='xr-wrap' style='display:none'><div class='xr-header'><div class='xr-obj-type'>xarray.DataArray</div><div class='xr-array-name'></div><ul class='xr-dim-list'><li><span class='xr-has-index'>FilterID</span>: 8162</li><li><span class='xr-has-index'>Te</span>: 82</li><li><span class='xr-has-index'>logg</span>: 13</li></ul></div><ul class='xr-sections'><li class='xr-section-item'><div class='xr-array-wrap'><input id='section-6111a7e5-74ab-4ad6-9089-e905f1cc3b79' class='xr-array-in' type='checkbox' checked><label for='section-6111a7e5-74ab-4ad6-9089-e905f1cc3b79' title='Show/hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-array-preview xr-preview'><span>...</span></div><div class='xr-array-data'><pre>[8700692 values with dtype=float64]</pre></div></div></li><li class='xr-section-item'><input id='section-f2377c58-a64b-4787-adaf-79383591a718' class='xr-section-summary-in' type='checkbox'  checked><label for='section-f2377c58-a64b-4787-adaf-79383591a718' class='xr-section-summary' >Coordinates: <span>(3)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><ul class='xr-var-list'><li class='xr-var-item'><div class='xr-var-name'><span class='xr-has-index'>FilterID</span></div><div class='xr-var-dims'>(FilterID)</div><div class='xr-var-dtype'>&lt;U38</div><div class='xr-var-preview xr-preview'>&#x27;IUE/IUE.1250-1300&#x27; ... &#x27;QUIJOTE...</div><input id='attrs-abbf4265-7f70-4c7d-aa21-33cd780d0708' class='xr-var-attrs-in' type='checkbox' disabled><label for='attrs-abbf4265-7f70-4c7d-aa21-33cd780d0708' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-5eee8901-357a-418a-8754-339bde241baa' class='xr-var-data-in' type='checkbox'><label for='data-5eee8901-357a-418a-8754-339bde241baa' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'></dl></div><div class='xr-var-data'><pre>array([&#x27;IUE/IUE.1250-1300&#x27;, &#x27;HST/STIS_FUV.F25LYA_G140M&#x27;,
           &#x27;HST/STIS_FUV.F25LYA_G140L&#x27;, ..., &#x27;QUIJOTE/MFI.13GHz_H1&#x27;,
           &#x27;QUIJOTE/MFI.11GHz_H1&#x27;, &#x27;QUIJOTE/MFI.11GHz_H3&#x27;], dtype=&#x27;&lt;U38&#x27;)</pre></div></li><li class='xr-var-item'><div class='xr-var-name'><span class='xr-has-index'>Te</span></div><div class='xr-var-dims'>(Te)</div><div class='xr-var-dtype'>int32</div><div class='xr-var-preview xr-preview'>5000 5250 5500 ... 70000 80000</div><input id='attrs-d86d195e-6fdf-4c86-bdb8-7e7c728b1ca3' class='xr-var-attrs-in' type='checkbox' disabled><label for='attrs-d86d195e-6fdf-4c86-bdb8-7e7c728b1ca3' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-9a745387-89b8-4754-9c17-d340c7da3724' class='xr-var-data-in' type='checkbox'><label for='data-9a745387-89b8-4754-9c17-d340c7da3724' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'></dl></div><div class='xr-var-data'><pre>array([ 5000,  5250,  5500,  5750,  6000,  6250,  6500,  6750,  7000,  7250,
            7500,  7750,  8000,  8250,  8500,  8750,  9000,  9250,  9500,  9750,
           10000, 10250, 10500, 10750, 11000, 11250, 11500, 11750, 12000, 12250,
           12500, 12750, 13000, 13250, 13500, 13750, 14000, 14250, 14500, 14750,
           15000, 15250, 15500, 15750, 16000, 16250, 16500, 16750, 17000, 17250,
           17500, 17750, 18000, 18250, 18500, 18750, 19000, 19250, 19500, 19750,
           20000, 21000, 22000, 23000, 24000, 25000, 26000, 27000, 28000, 29000,
           30000, 32000, 34000, 35000, 36000, 38000, 40000, 45000, 50000, 60000,
           70000, 80000], dtype=int32)</pre></div></li><li class='xr-var-item'><div class='xr-var-name'><span class='xr-has-index'>logg</span></div><div class='xr-var-dims'>(logg)</div><div class='xr-var-dtype'>float64</div><div class='xr-var-preview xr-preview'>6.5 6.75 7.0 7.25 ... 9.0 9.25 9.5</div><input id='attrs-5d7b45e3-6cf5-4882-9d87-f2eb72def909' class='xr-var-attrs-in' type='checkbox' disabled><label for='attrs-5d7b45e3-6cf5-4882-9d87-f2eb72def909' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-a67ad7a5-436d-4034-ba0d-d1323861d505' class='xr-var-data-in' type='checkbox'><label for='data-a67ad7a5-436d-4034-ba0d-d1323861d505' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'></dl></div><div class='xr-var-data'><pre>array([6.5 , 6.75, 7.  , 7.25, 7.5 , 7.75, 8.  , 8.25, 8.5 , 8.75, 9.  , 9.25,
           9.5 ])</pre></div></li></ul></div></li><li class='xr-section-item'><input id='section-da38729d-fe91-44f4-879a-8ca830d3a4c8' class='xr-section-summary-in' type='checkbox'  ><label for='section-da38729d-fe91-44f4-879a-8ca830d3a4c8' class='xr-section-summary' >Indexes: <span>(3)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><ul class='xr-var-list'><li class='xr-var-item'><div class='xr-index-name'><div>FilterID</div></div><div class='xr-index-preview'>PandasIndex</div><div></div><input id='index-ca5fa581-4386-4089-9f5b-0288c3a1f146' class='xr-index-data-in' type='checkbox'/><label for='index-ca5fa581-4386-4089-9f5b-0288c3a1f146' title='Show/Hide index repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-index-data'><pre>PandasIndex(Index([&#x27;IUE/IUE.1250-1300&#x27;, &#x27;HST/STIS_FUV.F25LYA_G140M&#x27;,
           &#x27;HST/STIS_FUV.F25LYA_G140L&#x27;, &#x27;TD1/SPEC.F1360&#x27;, &#x27;HST/STIS_FUV.F25LYA&#x27;,
           &#x27;TD1/SPEC.F1380&#x27;, &#x27;TD1/SPEC.F1400&#x27;, &#x27;TD1/TD1.140&#x27;, &#x27;TD1/NARROW1.1400&#x27;,
           &#x27;HST/ACS_SBC.F122M&#x27;,
           ...
           &#x27;COSMOSOMAS/COSMO15.17GHz&#x27;, &#x27;COSMOSOMAS/COSMO15.15GHz&#x27;,
           &#x27;QUIJOTE/MFI.17GHz_H2&#x27;, &#x27;COSMOSOMAS/COSMO15.13GHz&#x27;,
           &#x27;COSMOSOMAS/COSMO11.Ch1&#x27;, &#x27;COSMOSOMAS/COSMO11.Ch2&#x27;,
           &#x27;QUIJOTE/MFI.13GHz_H3&#x27;, &#x27;QUIJOTE/MFI.13GHz_H1&#x27;, &#x27;QUIJOTE/MFI.11GHz_H1&#x27;,
           &#x27;QUIJOTE/MFI.11GHz_H3&#x27;],
          dtype=&#x27;object&#x27;, name=&#x27;FilterID&#x27;, length=8162))</pre></div></li><li class='xr-var-item'><div class='xr-index-name'><div>Te</div></div><div class='xr-index-preview'>PandasIndex</div><div></div><input id='index-81225a49-8f62-4741-97d3-114e3ac73c91' class='xr-index-data-in' type='checkbox'/><label for='index-81225a49-8f62-4741-97d3-114e3ac73c91' title='Show/Hide index repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-index-data'><pre>PandasIndex(Index([ 5000,  5250,  5500,  5750,  6000,  6250,  6500,  6750,  7000,  7250,
            7500,  7750,  8000,  8250,  8500,  8750,  9000,  9250,  9500,  9750,
           10000, 10250, 10500, 10750, 11000, 11250, 11500, 11750, 12000, 12250,
           12500, 12750, 13000, 13250, 13500, 13750, 14000, 14250, 14500, 14750,
           15000, 15250, 15500, 15750, 16000, 16250, 16500, 16750, 17000, 17250,
           17500, 17750, 18000, 18250, 18500, 18750, 19000, 19250, 19500, 19750,
           20000, 21000, 22000, 23000, 24000, 25000, 26000, 27000, 28000, 29000,
           30000, 32000, 34000, 35000, 36000, 38000, 40000, 45000, 50000, 60000,
           70000, 80000],
          dtype=&#x27;int32&#x27;, name=&#x27;Te&#x27;))</pre></div></li><li class='xr-var-item'><div class='xr-index-name'><div>logg</div></div><div class='xr-index-preview'>PandasIndex</div><div></div><input id='index-67925933-2a58-4ad3-9a27-3123704ea1d0' class='xr-index-data-in' type='checkbox'/><label for='index-67925933-2a58-4ad3-9a27-3123704ea1d0' title='Show/Hide index repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-index-data'><pre>PandasIndex(Index([6.5, 6.75, 7.0, 7.25, 7.5, 7.75, 8.0, 8.25, 8.5, 8.75, 9.0, 9.25, 9.5], dtype=&#x27;float64&#x27;, name=&#x27;logg&#x27;))</pre></div></li></ul></div></li><li class='xr-section-item'><input id='section-25978b29-a554-4ffb-a8cd-b48f65e9e40b' class='xr-section-summary-in' type='checkbox'  checked><label for='section-25978b29-a554-4ffb-a8cd-b48f65e9e40b' class='xr-section-summary' >Attributes: <span>(3)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><dl class='xr-attrs'><dt><span>Wavelengths :</span></dt><dd>[1.28470318e+03 1.35336435e+03 1.35623969e+03 ... 3.86833493e+08
     5.62983439e+08 5.99076004e+08]</dd><dt><span>unit :</span></dt><dd>erg/s/cm2/A</dd><dt><span>long_name :</span></dt><dd>Flux</dd></dl></div></li></ul></div></div>



Creating VOSA input files
-------------------------

VOSA input format
~~~~~~~~~~~~~~~~~

VOSA requires following format for uploading the photometric information
files

+------+------+------+------+------+------+------+------+------+------+
| ob   | RA   | DEC  | dis  | Av   | fi   | flux | e    | pnt  | obj  |
| ject |      |      |      |      | lter |      | rror | opts | opts |
+======+======+======+======+======+======+======+======+======+======+
| —    | —    | —    | —    | —    | —    | —    | —    | —    | —    |
+------+------+------+------+------+------+------+------+------+------+

Identify the filters in the files you are uploading (create the
“filter_list” accordingly)

Create a file with “name, ra, dec, magnitudes and magnitude errors”.
This “photomety_file” be converted to VOSA format

-  Note: This code is for a cluster, hence distance and extinction is
   kept constant

The VOSA_input.txt file has magnitudes in the “flux” column. So while
uploading to VOSA, keep “file type” as “magnitudes”

.. code:: ipython3

    import pandas as pd
    import numpy as np
    
    df_photometry = pd.read_csv('data/example_photomety_file.csv', engine='python')
    df_photometry.set_index('name', inplace=True)
    
    # Adding distance and Av values
    distance      = '831.76+-11'       # [pc]
    Av            = '0.1736+-0.017'    # 3.1 * E(B-V) [mag]
    df_photometry['distance'] = distance
    df_photometry['Av'] = Av
    
    print((df_photometry.columns))
    df_photometry.head()


.. parsed-literal::

    Index(['ra', 'dec', 'B', 'e_B', 'V', 'e_V', 'I', 'e_I', 'U', 'e_U', 'R', 'e_R',
           'F148W', 'e_F148W', 'F154W', 'e_F154W', 'F169M', 'e_F169M', 'distance',
           'Av'],
          dtype='object')
    



.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>ra</th>
          <th>dec</th>
          <th>B</th>
          <th>e_B</th>
          <th>V</th>
          <th>e_V</th>
          <th>I</th>
          <th>e_I</th>
          <th>U</th>
          <th>e_U</th>
          <th>R</th>
          <th>e_R</th>
          <th>F148W</th>
          <th>e_F148W</th>
          <th>F154W</th>
          <th>e_F154W</th>
          <th>F169M</th>
          <th>e_F169M</th>
          <th>distance</th>
          <th>Av</th>
        </tr>
        <tr>
          <th>name</th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>WOCS2002</th>
          <td>132.8492</td>
          <td>11.83040</td>
          <td>12.300</td>
          <td>NaN</td>
          <td>11.540</td>
          <td>NaN</td>
          <td>11.060</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>18.66000</td>
          <td>0.058</td>
          <td>18.69300</td>
          <td>0.065</td>
          <td>18.58200</td>
          <td>0.058</td>
          <td>831.76+-11</td>
          <td>0.1736+-0.017</td>
        </tr>
        <tr>
          <th>WOCS1007</th>
          <td>132.8931</td>
          <td>11.85297</td>
          <td>11.070</td>
          <td>NaN</td>
          <td>10.960</td>
          <td>NaN</td>
          <td>11.072</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>17.03341</td>
          <td>0.018</td>
          <td>16.87171</td>
          <td>0.016</td>
          <td>16.63688</td>
          <td>0.019</td>
          <td>831.76+-11</td>
          <td>0.1736+-0.017</td>
        </tr>
        <tr>
          <th>Y1168</th>
          <td>132.8331</td>
          <td>11.81147</td>
          <td>18.409</td>
          <td>NaN</td>
          <td>18.632</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>18.393</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>17.02925</td>
          <td>0.039</td>
          <td>17.06582</td>
          <td>0.033</td>
          <td>17.17563</td>
          <td>0.033</td>
          <td>831.76+-11</td>
          <td>0.1736+-0.017</td>
        </tr>
      </tbody>
    </table>
    </div>



Renaming the columns with SVO Filter Profile Service names
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  http://svo2.cab.inta-csic.es/theory/fps/index.php?mode=browse

.. code:: ipython3

    filter_list = ['KPNO/Mosaic.B', 'KPNO/Mosaic.V','KPNO/Mosaic.I', 'KPNO/Mosaic.U','KPNO/Mosaic.R',
                   'Astrosat/UVIT.F148W','Astrosat/UVIT.F154W', 'Astrosat/UVIT.F169M']
    
    df_photometry.rename(columns={'ra'      : 'ra', 
                                  'dec'     : 'dec',
                                  'B'       : 'KPNO/Mosaic.B',
                                  'e_B'     : 'e_KPNO/Mosaic.B',
                                  'V'       : 'KPNO/Mosaic.V',
                                  'e_V'     : 'e_KPNO/Mosaic.V',
                                  'I'       : 'KPNO/Mosaic.I',
                                  'e_I'     : 'e_KPNO/Mosaic.I',
                                  'U'       : 'KPNO/Mosaic.U',
                                  'e_U'     : 'e_KPNO/Mosaic.U',
                                  'R'       : 'KPNO/Mosaic.R',
                                  'e_R'     : 'e_KPNO/Mosaic.R',
                                  'F148W'   : 'Astrosat/UVIT.F148W',
                                  'e_F148W' : 'e_Astrosat/UVIT.F148W',
                                  'F154W'   : 'Astrosat/UVIT.F154W',
                                  'e_F154W' : 'e_Astrosat/UVIT.F154W',
                                  'F169M'   : 'Astrosat/UVIT.F169M',
                                  'e_F169M' : 'e_Astrosat/UVIT.F169M',
                                  'distance': 'distance',
                                  'Av'      : 'Av'}, 
                         inplace=True)
    df_photometry.head()




.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>ra</th>
          <th>dec</th>
          <th>KPNO/Mosaic.B</th>
          <th>e_KPNO/Mosaic.B</th>
          <th>KPNO/Mosaic.V</th>
          <th>e_KPNO/Mosaic.V</th>
          <th>KPNO/Mosaic.I</th>
          <th>e_KPNO/Mosaic.I</th>
          <th>KPNO/Mosaic.U</th>
          <th>e_KPNO/Mosaic.U</th>
          <th>KPNO/Mosaic.R</th>
          <th>e_KPNO/Mosaic.R</th>
          <th>Astrosat/UVIT.F148W</th>
          <th>e_Astrosat/UVIT.F148W</th>
          <th>Astrosat/UVIT.F154W</th>
          <th>e_Astrosat/UVIT.F154W</th>
          <th>Astrosat/UVIT.F169M</th>
          <th>e_Astrosat/UVIT.F169M</th>
          <th>distance</th>
          <th>Av</th>
        </tr>
        <tr>
          <th>name</th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>WOCS2002</th>
          <td>132.8492</td>
          <td>11.83040</td>
          <td>12.300</td>
          <td>NaN</td>
          <td>11.540</td>
          <td>NaN</td>
          <td>11.060</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>18.66000</td>
          <td>0.058</td>
          <td>18.69300</td>
          <td>0.065</td>
          <td>18.58200</td>
          <td>0.058</td>
          <td>831.76+-11</td>
          <td>0.1736+-0.017</td>
        </tr>
        <tr>
          <th>WOCS1007</th>
          <td>132.8931</td>
          <td>11.85297</td>
          <td>11.070</td>
          <td>NaN</td>
          <td>10.960</td>
          <td>NaN</td>
          <td>11.072</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>17.03341</td>
          <td>0.018</td>
          <td>16.87171</td>
          <td>0.016</td>
          <td>16.63688</td>
          <td>0.019</td>
          <td>831.76+-11</td>
          <td>0.1736+-0.017</td>
        </tr>
        <tr>
          <th>Y1168</th>
          <td>132.8331</td>
          <td>11.81147</td>
          <td>18.409</td>
          <td>NaN</td>
          <td>18.632</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>18.393</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>NaN</td>
          <td>17.02925</td>
          <td>0.039</td>
          <td>17.06582</td>
          <td>0.033</td>
          <td>17.17563</td>
          <td>0.033</td>
          <td>831.76+-11</td>
          <td>0.1736+-0.017</td>
        </tr>
      </tbody>
    </table>
    </div>



Creating the VOSA compliant input file from the dataframe
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    def create_VOSA_input_file(df_photometry, filter_list, output_name=None):
        # combining data from all stars to make the VOSA upload file 
        df_VOSA_input = pd.DataFrame(columns = ['object', 'ra', 'dec','distance','Av','filter','flux','error','pntopts','objopts']) 
        counter=0
        for name in df_photometry.index:
            for filter_name in filter_list:
                df_VOSA_input.loc[counter, 'object']   = name
                df_VOSA_input.loc[counter, 'ra']       = df_photometry['ra'][name]
                df_VOSA_input.loc[counter, 'dec']      = df_photometry['dec'][name]
                df_VOSA_input.loc[counter, 'distance'] = df_photometry['distance'][name]
                df_VOSA_input.loc[counter, 'Av']       = df_photometry['Av'][name]
                df_VOSA_input.loc[counter, 'filter']   = filter_name
                df_VOSA_input.loc[counter, 'flux']     = df_photometry[filter_name][name]
                df_VOSA_input.loc[counter, 'error']    = df_photometry['e_'+filter_name][name]
                df_VOSA_input.loc[counter, 'pntopts']  = '---'
                df_VOSA_input.loc[counter, 'objopts']  = '---'
    
                # Rewriting the row if flux is not available for this star-filter combination
                if ~np.isnan(df_photometry[filter_name][name]):
                    counter=counter+1
        df_VOSA_input.fillna('---', inplace=True)
        if output_name !=None:
            df_VOSA_input.to_csv(output_name, header=None, index=None, sep=' ')
        return df_VOSA_input
    
    create_VOSA_input_file(df_photometry, filter_list, output_name='data/example_VOSA_input_file.txt')


.. parsed-literal::

    /var/folders/2v/ztlgg9951mdc3j88_z3gbscc0000gn/T/ipykernel_5393/3182199717.py:21: FutureWarning: Downcasting object dtype arrays on .fillna, .ffill, .bfill is deprecated and will change in a future version. Call result.infer_objects(copy=False) instead. To opt-in to the future behavior, set `pd.set_option('future.no_silent_downcasting', True)`
      df_VOSA_input.fillna('---', inplace=True)
    



.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>object</th>
          <th>ra</th>
          <th>dec</th>
          <th>distance</th>
          <th>Av</th>
          <th>filter</th>
          <th>flux</th>
          <th>error</th>
          <th>pntopts</th>
          <th>objopts</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>0</th>
          <td>WOCS2002</td>
          <td>132.8492</td>
          <td>11.83040</td>
          <td>831.76+-11</td>
          <td>0.1736+-0.017</td>
          <td>KPNO/Mosaic.B</td>
          <td>12.30000</td>
          <td>---</td>
          <td>---</td>
          <td>---</td>
        </tr>
        <tr>
          <th>1</th>
          <td>WOCS2002</td>
          <td>132.8492</td>
          <td>11.83040</td>
          <td>831.76+-11</td>
          <td>0.1736+-0.017</td>
          <td>KPNO/Mosaic.V</td>
          <td>11.54000</td>
          <td>---</td>
          <td>---</td>
          <td>---</td>
        </tr>
        <tr>
          <th>2</th>
          <td>WOCS2002</td>
          <td>132.8492</td>
          <td>11.83040</td>
          <td>831.76+-11</td>
          <td>0.1736+-0.017</td>
          <td>KPNO/Mosaic.I</td>
          <td>11.06000</td>
          <td>---</td>
          <td>---</td>
          <td>---</td>
        </tr>
        <tr>
          <th>3</th>
          <td>WOCS2002</td>
          <td>132.8492</td>
          <td>11.83040</td>
          <td>831.76+-11</td>
          <td>0.1736+-0.017</td>
          <td>Astrosat/UVIT.F148W</td>
          <td>18.66000</td>
          <td>0.058</td>
          <td>---</td>
          <td>---</td>
        </tr>
        <tr>
          <th>4</th>
          <td>WOCS2002</td>
          <td>132.8492</td>
          <td>11.83040</td>
          <td>831.76+-11</td>
          <td>0.1736+-0.017</td>
          <td>Astrosat/UVIT.F154W</td>
          <td>18.69300</td>
          <td>0.065</td>
          <td>---</td>
          <td>---</td>
        </tr>
        <tr>
          <th>5</th>
          <td>WOCS2002</td>
          <td>132.8492</td>
          <td>11.83040</td>
          <td>831.76+-11</td>
          <td>0.1736+-0.017</td>
          <td>Astrosat/UVIT.F169M</td>
          <td>18.58200</td>
          <td>0.058</td>
          <td>---</td>
          <td>---</td>
        </tr>
        <tr>
          <th>6</th>
          <td>WOCS1007</td>
          <td>132.8931</td>
          <td>11.85297</td>
          <td>831.76+-11</td>
          <td>0.1736+-0.017</td>
          <td>KPNO/Mosaic.B</td>
          <td>11.07000</td>
          <td>---</td>
          <td>---</td>
          <td>---</td>
        </tr>
        <tr>
          <th>7</th>
          <td>WOCS1007</td>
          <td>132.8931</td>
          <td>11.85297</td>
          <td>831.76+-11</td>
          <td>0.1736+-0.017</td>
          <td>KPNO/Mosaic.V</td>
          <td>10.96000</td>
          <td>---</td>
          <td>---</td>
          <td>---</td>
        </tr>
        <tr>
          <th>8</th>
          <td>WOCS1007</td>
          <td>132.8931</td>
          <td>11.85297</td>
          <td>831.76+-11</td>
          <td>0.1736+-0.017</td>
          <td>KPNO/Mosaic.I</td>
          <td>11.07200</td>
          <td>---</td>
          <td>---</td>
          <td>---</td>
        </tr>
        <tr>
          <th>9</th>
          <td>WOCS1007</td>
          <td>132.8931</td>
          <td>11.85297</td>
          <td>831.76+-11</td>
          <td>0.1736+-0.017</td>
          <td>Astrosat/UVIT.F148W</td>
          <td>17.03341</td>
          <td>0.018</td>
          <td>---</td>
          <td>---</td>
        </tr>
        <tr>
          <th>10</th>
          <td>WOCS1007</td>
          <td>132.8931</td>
          <td>11.85297</td>
          <td>831.76+-11</td>
          <td>0.1736+-0.017</td>
          <td>Astrosat/UVIT.F154W</td>
          <td>16.87171</td>
          <td>0.016</td>
          <td>---</td>
          <td>---</td>
        </tr>
        <tr>
          <th>11</th>
          <td>WOCS1007</td>
          <td>132.8931</td>
          <td>11.85297</td>
          <td>831.76+-11</td>
          <td>0.1736+-0.017</td>
          <td>Astrosat/UVIT.F169M</td>
          <td>16.63688</td>
          <td>0.019</td>
          <td>---</td>
          <td>---</td>
        </tr>
        <tr>
          <th>12</th>
          <td>Y1168</td>
          <td>132.8331</td>
          <td>11.81147</td>
          <td>831.76+-11</td>
          <td>0.1736+-0.017</td>
          <td>KPNO/Mosaic.B</td>
          <td>18.40900</td>
          <td>---</td>
          <td>---</td>
          <td>---</td>
        </tr>
        <tr>
          <th>13</th>
          <td>Y1168</td>
          <td>132.8331</td>
          <td>11.81147</td>
          <td>831.76+-11</td>
          <td>0.1736+-0.017</td>
          <td>KPNO/Mosaic.V</td>
          <td>18.63200</td>
          <td>---</td>
          <td>---</td>
          <td>---</td>
        </tr>
        <tr>
          <th>14</th>
          <td>Y1168</td>
          <td>132.8331</td>
          <td>11.81147</td>
          <td>831.76+-11</td>
          <td>0.1736+-0.017</td>
          <td>KPNO/Mosaic.U</td>
          <td>18.39300</td>
          <td>---</td>
          <td>---</td>
          <td>---</td>
        </tr>
        <tr>
          <th>15</th>
          <td>Y1168</td>
          <td>132.8331</td>
          <td>11.81147</td>
          <td>831.76+-11</td>
          <td>0.1736+-0.017</td>
          <td>Astrosat/UVIT.F148W</td>
          <td>17.02925</td>
          <td>0.039</td>
          <td>---</td>
          <td>---</td>
        </tr>
        <tr>
          <th>16</th>
          <td>Y1168</td>
          <td>132.8331</td>
          <td>11.81147</td>
          <td>831.76+-11</td>
          <td>0.1736+-0.017</td>
          <td>Astrosat/UVIT.F154W</td>
          <td>17.06582</td>
          <td>0.033</td>
          <td>---</td>
          <td>---</td>
        </tr>
        <tr>
          <th>17</th>
          <td>Y1168</td>
          <td>132.8331</td>
          <td>11.81147</td>
          <td>831.76+-11</td>
          <td>0.1736+-0.017</td>
          <td>Astrosat/UVIT.F169M</td>
          <td>17.17563</td>
          <td>0.033</td>
          <td>---</td>
          <td>---</td>
        </tr>
      </tbody>
    </table>
    </div>



Now upload the file at
http://svo2.cab.inta-csic.es/theory/vosa/index.php?action=myfiles&otype=star&seeall=

Make sure to change the File type to \`magnitude`\`

Keep SED Type: Flux vs Lambda

Select the file and search through VO for all possible detections Look
at the SEDs,

Possibly remove some telescopes (e.g. SDSS photometry is problematic for
bright M67 stars, Gaia DR3 synthetic photometry may or may not be
usable)

Create flux corrected files from VOSA
-------------------------------------

-  DEPRECATED

   -  NO LONGER REQUIRED AS ``bsf.load_data()`` CAN NOW USE VOSA FILES
      DIRECTLY

.. code:: ipython3

    import pandas as pd
    def create_extinction_corrected_files_from_VOSA(name, file_name, output_name, verbose=True):
        # Reading parameters derived by VOSA
        flux         = pd.read_csv(file_name, engine='python', comment='#', delim_whitespace= True, skipinitialspace=True, header=None)
        flux.columns = ['FilterID','wavelength','Obs.Flux','Obs.Error','flux','error','model_flux_A','Fitted','Excess','FitExc','UpLim']
    
        # Removing filters noted as "upper limit"
        if verbose: print('WARNING: %s removed due to upper limit'%(flux[flux['UpLim']=='1']['FilterID'].values))
        flux         = flux[flux['UpLim'] == '---']
        flux         = flux.drop(columns=['Obs.Flux','Obs.Error','Fitted','Excess','FitExc','UpLim','model_flux_A'])
        flux.to_csv(output_name, index=False)
        return flux
    
    DIR_OBS = 'data/'
    
    for idx, name in enumerate(['WOCS2002','WOCS1007','Y1168']):    
        print(name)
        # file_name = DIR_OBS + name +'/bestfitp/'+ name +'.bfit.phot.dat'
        file_name = 'data/vosa_results/%s.bfit.phot.dat'%name
        flux = create_extinction_corrected_files_from_VOSA(name,file_name, output_name='data/extinction_corrected_flux_files/'+name+'.csv')


.. parsed-literal::

    WOCS2002
    WARNING: ['WISE/WISE.W4'] removed due to upper limit
    WOCS1007
    WARNING: ['WISE/WISE.W4'] removed due to upper limit
    Y1168
    WARNING: ['WISE/WISE.W3' 'WISE/WISE.W4'] removed due to upper limit
    

.. parsed-literal::

    /var/folders/2v/ztlgg9951mdc3j88_z3gbscc0000gn/T/ipykernel_5393/1800319444.py:4: FutureWarning: The 'delim_whitespace' keyword in pd.read_csv is deprecated and will be removed in a future version. Use ``sep='\s+'`` instead
      flux         = pd.read_csv(file_name, engine='python', comment='#', delim_whitespace= True, skipinitialspace=True, header=None)
    /var/folders/2v/ztlgg9951mdc3j88_z3gbscc0000gn/T/ipykernel_5393/1800319444.py:4: FutureWarning: The 'delim_whitespace' keyword in pd.read_csv is deprecated and will be removed in a future version. Use ``sep='\s+'`` instead
      flux         = pd.read_csv(file_name, engine='python', comment='#', delim_whitespace= True, skipinitialspace=True, header=None)
    /var/folders/2v/ztlgg9951mdc3j88_z3gbscc0000gn/T/ipykernel_5393/1800319444.py:4: FutureWarning: The 'delim_whitespace' keyword in pd.read_csv is deprecated and will be removed in a future version. Use ``sep='\s+'`` instead
      flux         = pd.read_csv(file_name, engine='python', comment='#', delim_whitespace= True, skipinitialspace=True, header=None)
    

Logo
----

.. code:: ipython3

    import xarray as xr
    import src.binary_sed_fitting as bsf
    import matplotlib.pyplot as plt
    
    da = xr.open_dataarray(bsf.DIR_MODELS+'koester_spectral_library.nc')
    y1 = da.sel(Te=6000).sel(logg=7.0)*4e2
    y2 = da.sel(Te=15000).sel(logg=7.0)*1.5
    
    plt.style.use('dark_background')
    fig, ax = plt.subplots(figsize=(2,1))
    
    ax.fill_between(da.Wavelength, y1, alpha=0.6, color='#FF6961')
    ax.fill_between(da.Wavelength, y2, alpha=0.6, color='#6198ff')
    
    ax.loglog()
    ax.set_xlim(969,8000)
    ax.set_ylim(3000000)
    ax.axis('off')
    ax.text(0.5, 0.5, 'Binary     \n  SED      \n    Fitting', fontsize=14, transform=ax.transAxes,
            va='center', ha='center', color='white', family='serif', variant='small-caps', weight='bold')
    plt.savefig('plots/logo.png', dpi=200, bbox_inches='tight')



.. image:: examples_io/output_19_0.png

