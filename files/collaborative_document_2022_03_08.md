![](https://i.imgur.com/iywjz8s.png)


# 2022-03-07 Astronomical Data Science (day 2).

Welcome to The Workshop Collaborative Document.

This Document is synchronized as you type, so that everyone viewing this page sees the same text. This allows you to collaborate seamlessly on documents.

----------------------------------------------------------------------------

This is the Document for today: [link](https://hackmd.io/@O0tsDNPbTlyhyGiiCMaLIw/B1JW-KNZc)

Collaborative Document day 1: [link](https://hackmd.io/@O0tsDNPbTlyhyGiiCMaLIw/SkeAMVmZ9)

Collaborative Document day 2: [link](https://hackmd.io/@O0tsDNPbTlyhyGiiCMaLIw/B1JW-KNZc)

Collaborative Document day 3: [link](<url>)


## 👮Code of Conduct

* Participants are expected to follow those guidelines:
* Use welcoming and inclusive language.
* Be respectful of different viewpoints and experiences.
* Gracefully accept constructive criticism.
* Focus on what is best for the community.
* Show courtesy and respect towards other community members.
 
## ⚖️ License

All content is publicly available under the Creative Commons Attribution License: [creativecommons.org/licenses/by/4.0/](https://creativecommons.org/licenses/by/4.0/).

## 🙋Getting help

To ask a question, type `/hand` in the chat window.

To get help, type `/help` in the chat window.

You can ask questions in the document or chat window and helpers will try to help you.


```python=
phi2 = results_df['phi2']
phi2_min = -1.0 * u.degree
phi2_max = 1.0 * u.degree
phi2 < phi2_max
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

## 🖥 Workshop website

[link](<url>)

🛠 Setup

[link](<url>)

Download files

[link](<url>)

## 👩‍🏫👩‍💻🎓 Instructors
Hanno Spreeuw, Johan Hidding

## 🧑‍🙋 Helpers
Dafne van Kuppevelt, Suvayu Ali

## 👩‍💻👩‍💼👨‍🔬🧑‍🔬🧑‍🚀🧙‍♂️🔧 Roll Call
Name/ pronouns (optional) / job, role / social media (twitter, github, ...) / background or interests (optional) / cit



## 🗓️ Agenda
|time|what|
|---|---|
| 9:00 | Welcome and Icebreaker |
| 9:15 | Plotting and Pandas |
|10:00 | Transform and Select |
|10:15 | break |
|10:30 | Transform and Select, cont'd |
|11:30 | break |
|11:45 | JOIN |
|12:45 | Closing and feedback |

:
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

## 🧠 Collaborative Notes
Recap of last part (plotting) of yesterday, see [notes of yesterday](https://hackmd.io/BGoihtChQgazQjA-9S2lEQ?view#Plotting).

Let's transform to the GD-1 frame:
```python=
transformed = skycoord.transform_to(gd1_frame)
```

We are going to do a strict selection on position in sky, and proper motion.
We need to do reflex correction, to correct for the motion of our solar system around the galactic center. The algorithm for this is already available in astropy (see [documentation](https://gala-astro.readthedocs.io/en/latest/api/gala.coordinates.reflex_correct.html)).
```python=
from gala.coordinates import reflex_correct

skycoord_gd1 = reflex_correct(transformed)
```

Let's see if we indeed have a rectangle in GD-1.
```python=
x = skycoord_gd1.phi1
y = skycoord_gd1.phi2
plt.plot(x, y, 'ko', markersize=0.1, alpha=0.1)

plt.xlabel('phi1 (degree GD1)')
plt.ylabel('phi2 (degree GD1)')
```

We are going to add columns to our dataset, which will help to select the right data. Note that polygon_results currently 'lives' in ICRS with `ra` and `dec`. However, we are going to add columns `phi1` and `phi2` that 'live' in GD-1. A bit dirty, but it's going to be useful.

For that purpose, we will convert our data to a Pandas dataframe. Pandas is more generic than astropy. But we will lose metadata (e.g. units).
```python=
type(polygon_results)  # This is an astropy table

type(skycoord) # This is an astropy SkyCoord object

polygon_results['phi1'] = skycoord_gd1.phi1
polygon_results['phi2'] = skycoord_gd1.phi2
polygon_results['pm_phi1'] = skycoord_gd1.pm_phi1_cosphi2
polygon_results['pm_phi2'] = skycoord_gd1.pm_phi2

# Convert to a pandas dataframe
import pandas as pd
results_df = polygon_results.to_pandas()
results_df.shape
results_df.head() # Shows first five lines

# Create the function (as in the exercise)
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
Now use the function, and save the result. Because we use pandas, we can save as HDF5 and add tables to the same file later.
```python=
results_df = make_dataframe(polygon_results)
filename = 'gd1_data.hdf'
# `results_df` is the key for this table in the file:
results_df.to_hdf(filename, 'results_df', 'w')  # 'w': overwrite if already exists
```

Now we are going to use the table to recreate the plots from the papers. First we explore the data a bit:
```python=
results.describe() # outputs a table with column statistics: count, mean, quartile, etc

# Or if we want to see the units:
polygon_results.info()
```
Note that parallax is sometimes negative, which is not physically possible. 
These are measurement errors. We don't use the parallex
That's why we used the constant distance from reflex correction.

Now we plot the proper motion
```python=
x = results_df['pm_phi1']
y = results_df['pm_phi2']
plt.plot(x, y, 'ko', markersize=0.1, alpha=0.1)
    
plt.xlabel('Proper motion phi1 (mas/yr GD1 frame)')
plt.ylabel('Proper motion phi2 (mas/yr GD1 frame)'
```

We zoom in a bit more:
```python=
x = results_df['pm_phi1']
y = results_df['pm_phi2']
plt.plot(x, y, 'ko', markersize=0.1, alpha=0.3)
    
plt.xlabel('Proper motion phi1 (mas/yr GD1 frame)')
plt.ylabel('Proper motion phi2 (mas/yr GD1 frame)')

plt.xlim(-10, 10)
plt.ylim(-10, 10)
```

We want to get rid of the foreground stars, so we limit the window in `phi2`.
```python=
phi2_min = -1.0 * u.degree
phi2_max = 1.0 * u.degree

# Create a mask to be applied to phi2:
phi2 = results_df['phi2']
mask = (phi2 > phi2_min) & (phi2 < phi2_max)

# To see how many stars are selected:
mask.sum()
```

Let's put the selected stars in a pandas dataframe

```python=
centerline_df = results_df[mask]
len(centerline_df) # should match the count from above
len(centerline_df)/len(results_df) # fraction of selected stars

pm1_min = -8.9
pm1_max = -6.9
pm2_min = -2.2
pm2_max = 1.0

pm1_rect, pm2_rect = make_rectangle(pm1_min, pm1_max, pm2_min, pm2_max)
plot_proper_motion(centerline_df)
plt.plot(pm1_rect, pm2_rect, "-") # plot
```

Define a between function for further selection:
```python=
def between(series, low, high):
    return(series > low) & (series < high)

pm1 = results_df["pm_phi1"]
pm2 = results_df["pm_phi2"]

pm_mask = (between(pm1, pm1_min, pm1_max) & between(pm2, pm2_min, pm2_max))

pm_mask.sum() # returns the count after this second selection (~1000)

selected_df = results_df[pm_mask]
```
Note: using masks on the dataset as we are doing, this is not supported in `astropy`; another argument for using `pandas` for some of the analysis

Plot the selected stars and see if the stream can be identified:
```python=
def plot_pm_selection(df):
    "Plot in GD-1 spatial coordinates the location of the stars selected by proper motion"
    
    x = df["phi1"]
    y = df["phi2"]
    
    plt.plot(x, y, "ko", markersize=0.3, alpha=0.3)
    plt..xlabel("phi1 [deg]")
    plt..ylabel("phi2 [deg]")
    plt.title("Proper motion selection", fontsize="medium")
    
    plt.axis("equal")
    
plot_pm_selection(selected_df)
```

To write to the HDF5 file, do as before but be careful not to overwrite the file; use the append mode (`a` -> append, `w` -> overwrite).  Append is the pandas default, so no need to be too scared of overwriting your existing work.

```python=
selected_df.to_hdf(filename, "selected_df")
centerline_df.to_hdf(filename, "centerline_df")
```

How to read the HDF5 files and check what you have written:

```python=
# Python context manager, recommended for working with files
with pd.HDFStore(filename) as hdf: 
    print(pdf.keys())  # prints out all the keys stored in the file
    hdf.info()  # prints the stored dataframe keys, and their shape
```

To get a estimate of the file size:
```python=
from os.path import getsize

MB = 1024*1024
getsize(filename)/MB  # approximately 15 MB in this case
```

Can we do the same selection by writing a query and submitting to the Gaia database? The query has to be in right ascension/declination coordinates, it won't be a rectangle, but a more complicated polygon.

```python=
plot_proper_motion(centerline_df)

plt.plot(pm1_rect, pm2_rect)

x = selected_df["pm_phi1"]
y = selected_df["pm_phi2"]
plt.plot(x, y, "gx", markersize=0.3, alpha=0.3)
```

```python=
x = centerline_df["pmra"]
y = centerline_df["pmdec"]

plt.plot(x, y, "ko", ms=0.3, alpha=0.3)

x = selected_df["pmra"]
y = selected_df["pmdec"]

plt.plot(x, y, "gx", ms=0.3, alpha=0.3)

plt.xlabel("proper motion along right ascension (mas/yr)")
plt.ylabel("proper motion along declination (mas/yr)")

plt.xlim([-10, 4])
plt.ylim([-15, 2.5])
```

Note that the green rectangle has now transformed into an irregular shape.  This has to be covered if we are to send a query to the Gaia database

Let us use a `scipy` routine to make this easier (convex hull).
```python=
import numpy as np

# b/c scipy routine works on numpy arrays
points = selected_df[["pmra", "pmdec"]].to_numpy() 
# note: for older versions of pandas, you might need to do:
# dataframe.values instead

points.shape # (1049, 2)
```

```python=
from scipy.spatial import ConvexHull

hull = ConvexHull(points)
hull

hull.vertices # indices of the stars that are at the boundary (convex hull)
```

These are just the indices to the stars in the final selection, we need to convert it to proper motion coordinates in right ascension and declination.

```python=
pm_vertices = points[hull.vertices]
pm_vertices
```

Let us separate the coordinates into separate arrays:
```python=
pmra_poly, pmdec_poly = np.transpose(pm_vertices)
```
The above line transposes the above array from before: `N x 2` -> `2 x N`, `N` is the number of stars on the boundary

```python=
x = centerline_df['pmra']
y = centerline_df['pmdec']
plt.plot(x, y, 'ko', markersize=0.3, alpha=0.3)

x = selected_df['pmra']
y = selected_df['pmdec']
plt.plot(x, y, 'gx', markersize=0.3, alpha=0.3)

plt.plot(pmra_poly, pmdec_poly)
    
plt.xlabel('Proper motion phi1 (ICRS frame)')
plt.ylabel('Proper motion phi2 (ICRS frame)')

plt.xlim([-10, 0])
plt.ylim([-16, -10]);
```

Define the boundaries and transform:
```python=
phi1_min = -70 * u.degree
phi1_max = -20 * u.degree
phi2_min = -5 * u.degree
phi2_max = 5 * u.degree

phi1_rect, phi2_rect = make_rectangle(phi1_min, phi1_max, phi2_min, phi2_max)

corners = SkyCoord(phi1=phi1_rect, phi2=phi2_rect, frame=gd1_frame)
corners_icrs = corners.transform_to("icrs") # coordinate transform
```
Recall from yesterday, you need something more. This needs to be converted to string to be able to query

```python=
def skycoord_to_string(skycord):  # from day 1
    """Convert a 1-d list of SkyCoord to string in ADQL format"""
    return ' '.join(skycord.to_string()).replace(' ', ',')

sky_point_list = skycoord_to_string(corners_icrs)
# this should look something like: 135.306,8.39862,126.51,...
# 10 numbers, 4 corners, with the first point repeated at the end
```

Query:
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
s = np.array2string(pm_vertices.flatten(), max_line_width=1000, separator=",")
# same as above, list of the coordinates concatenated into a string
# note this includes square brackets at the ends, they need to be stripped
pm_point_list = s.strip("[]")  # ADQL should be okay with this
```

Updated query:
```python=
candidate_coord_query_base = """SELECT
{columns}
FROM gaiadr2.gaia_source
WHERE parallax < 1
  AND bp_rp BETWEEN -0.75 AND 2 
  AND 1 = CONTAINS(POINT(ra, dec), 
                   POLYGON({sky_point_list}))
  AND 1 = CONTAINS(POINT(pmra, pmdec), 
                   POLYGON({pm_point_list}))
"""

candidate_coord_pm_query = candidate_coord_query_base.format(
    columns=columns,
    sky_point_list=sky_point_list,
    pm_point_list=pm_point_list
)
```
Print out the query to get a feel how it looks:
```python
print(candidate_coord_pm_query)
```

You can now launch the query on the Gaia database. Make sure it is an asynchronous job, as it's a large job.
```python
candidate_coord_pm_job = Gaia.launch_job_async(candidate_coord_pm_query)
print(candidate_coord_pm_job)
```

Note: if you get an `HTTPError` (`Error 500`), this is probably an issue at the Gaia database end, maybe you can workaround the issue by using a circular boundary instead of a polygon.

## 📚 Resources

## Tips
(feedback to improve future iterations)

- For me the pace of today could have been a little bit quicker.

## Tops
(what you liked in today's session)

Nice to keep on plotting continuously, to check what we are actually doing/selecting

Nice course.
The idea that we can get only the selected data by refining the search and doing the computation on server itself was quite nice. 

Learning more complicated ADQL queries was nice and very usefull for the future.

Love those plots :)

I like being able to follow along with the datacarpentry site
