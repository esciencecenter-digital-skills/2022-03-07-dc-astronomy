# Astronomical Data Science Carpentry exercise thingy

List of participants for copy pasting:

List of break-out rooms:
* Room 1:
* Room 2:
* Room 3:
* Room 4:

```python=

```

## 🔧 Exercises
To check: do you have the object `skycoord` ready?
Also make sure to copy this code:
```python=
def make_rectangle(x1, x2, y1, y2):
    """Return the corners of a rectangle."""
    xs = [x1, x1, x2, x2, x1] # coordinates in phi1
    ys = [y1, y2, y2, y1, y1] # coordinates in phi2
    return xs, ys
```


### Exercise: plot in the new frame
Make a plot, similar as yesterdays plot, in the GD-1 frame.


Solution:
```python=
x = skycoord_gd1.phi1
y = skycoord_gd1.phi2
plt.plot(x, y, 'ko', markersize=0.1, alpha=0.1)

plt.xlabel('phi1 (degree GD1)')
plt.ylabel('phi2 (degree GD1)')
```

### Exercise make_dataframe
Create a function `make_dataframe` that takes as an input an astropy table, transforms to GD1, and outputs a pandas dataframe with the extra columns.
You may take the distance and radial_velocity as hard coded constants.


Solution:
```python=
def make_dataframe(table):
    """Transform coordinates from ICRS to GD-1 frame.
    
    table: Astropy Table
    
    returns: Pandas DataFrame
    """
    #Create a SkyCoord object with the coordinates and proper motions
    # in the input table
    skycoord = SkyCoord(
               ra=table['ra'], 
               dec=table['dec'],
               pm_ra_cosdec=table['pmra'],
               pm_dec=table['pmdec'], 
               distance=8*u.kpc, 
               radial_velocity=0*u.km/u.s)

    # Define the GD-1 reference frame
    gd1_frame = GD1Koposov10()

    # Transform input coordinates to the GD-1 reference frame
    transformed = skycoord.transform_to(gd1_frame)

    # Correct GD-1 coordinates for solar system motion around galactic center
    skycoord_gd1 = reflex_correct(transformed)

    #Add GD-1 reference frame columns for coordinates and proper motions
    table['phi1'] = skycoord_gd1.phi1
    table['phi2'] = skycoord_gd1.phi2
    table['pm_phi1'] = skycoord_gd1.pm_phi1_cosphi2
    table['pm_phi2'] = skycoord_gd1.pm_phi2

    # Create DataFrame
    df = table.to_pandas()

    return df
```

## Exercise: zoom in
Take this plot of the proper motion:
```python=
x = results_df['pm_phi1']
y = results_df['pm_phi2']
plt.plot(x, y, 'ko', markersize=0.1, alpha=0.1)
    
plt.xlabel('Proper motion phi1 (mas/yr GD1 frame)')
plt.ylabel('Proper motion phi2 (mas/yr GD1 frame)')
```
But zoom in a bit more on the 'blob'


## Exercise: plot properline motion of `centerline_df`
Since there are fewer stars, maybe you won't need to zoom as much as before
```python=
phi2_min = -1.0 * u.degree
phi2_max = 1.0 * u.degree

# Create a mask to be applied to phi2:
phi2 = results_df['phi2']
mask = (phi2 > phi2_min) & (phi2 < phi2_max)

centerline_df = results_df[mask]
```

```python=
x = centerline_pd['pm_phi1']
y = centerline_pd['pm_phi2']
plt.plot(x,y, 'ko', markersize=0.1, alpha=0.3)
plt.xlim(-10, 10)
plt.ylim(-10, 10)
plt.xlabel('motion phi1 (mas/yr GD1 frame)')
plt.ylabel('Proper motion phi2 (mas/yr GD1 frame)')
```

## Exercise: plot the `selected_df`
Plot the 1049 stars in the `selected_df` and see if you can see a band of stars in the sky
```python=
# `between` has been defined below in the notes
pm_mask = (between(pm1, pm1_min, pm1_max) & between(pm2, pm2_min, pm2_max))
selected_df = results_df[pm_mask]
```

```python=
x1 = selected_pd['phi1']
y1 = selected_pd['phi2']
plt.plot(x1,y1, 'ko', markersize=1, alpha=0.5)
plt.xlabel('phi1 (mas/yr GD1 frame)')
plt.ylabel('phi2 (mas/yr GD1 frame)')
```

## Exercise: add the two new dataframes you have now to the HDF5 file
Note: do not overwrite your existing file
```python=
dataframe.to_hdf(..., mode="w")  # overwrite
dataframe.to_hdf(..., mode="a")  # append
```

```python=
filename = 'gd1_data.hdf'
selected_df.to_hdf(filename, 'selected_df', mode='a')
centerline_df.to_hdf(filename, 'centerline_df', mode='a')
```


## Exercise: plot the figure using `pmra` & `pmdec` coordinates
```python=
x = centerline_df["pmra"]
y = centerline_df["pmdec"]

# plot commands go here

x = selected_df["pmra"]
y = selected_df["pmdec"]

# plot commands go here
```


## Exercise: check if the convex hull routine actually do the job?
Draw a green line along the boundary, and check by drawing a plot

Solution:
```python=
x = centerline_df['pmra']
y = centerline_df['pmdec']
plt.plot(x, y, 'ko', markersize=0.5, alpha=0.3)

x = selected_df['pmra']
y = selected_df['pmdec']
plt.plot(x, y, 'go', markersize=0.6, alpha=0.4)

plt.plot(pmra_poly, pmdec_poly, 'g-', alpha=0.7)

plt.xlabel('Proper motion RA (mas/yr ICRS frame)')
plt.ylabel('Proper motion DEC (mas/yr ICRS frame)')
plt.xlim(-10, 4)
plt.ylim(-16, 2)
```
![](https://i.imgur.com/zZE620H.png)

## Exercise: expand the base `candidate_coord_query` with an entry for the `pm_point_list`
```python=
columns = "source_id, ra, dec, pmra, pmdec"

candidate_coord_query_base = """SELECT
{columns}
FROM gaiadr2.gaia_source
WHERE parallax < 1
  AND bp_rp BETWEEN -0.75 AND 2 
  AND 1 = CONTAINS(POINT(ra, dec), 
                   POLYGON({sky_point_list}))
"""

candidate_coord_query = candidate_coord_query_base.format(
    columns=columns,
    sky_point_list=sky_point_list
)
```


```python=
candidate_coord_pm_query_base = """SELECT
{columns}
FROM gaiadr2.gaia_source
WHERE parallax < 1 
AND bp_rp BETWEEN -0.75 AND 2
AND 1 = CONTAINS(POINT(ra,dec),
                 POLYGON({sky_point_list}))
AND 1 = CONTAINS(POINT(pmra,pmdec),
                 POLYGON({pm_point_list}))"""
```
