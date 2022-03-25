![](https://i.imgur.com/iywjz8s.png)


# 2022-03-07 Astronomical Data Science (Day 1)

Welcome to The Workshop Collaborative Document.

This Document is synchronized as you type, so that everyone viewing this page sees the same text. This allows you to collaborate seamlessly on documents.

----------------------------------------------------------------------------

This is the Document for today: [link](https://hackmd.io/@O0tsDNPbTlyhyGiiCMaLIw/SkeAMVmZ9)

Collaborative Document day 1: [link](https://hackmd.io/@O0tsDNPbTlyhyGiiCMaLIw/SkeAMVmZ9)

Collaborative Document day 2: [link](https://hackmd.io/t25qWXHRRqWgvN6SVwcvAw)

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

## 🖥 Workshop website

[link](https://esciencecenter-digital-skills.github.io/2022-03-07-dc-astronomy/)

🛠 Setup

[link](https://datacarpentry.org/astronomy-python/setup.html)

Download files

[link](https://zenodo.org/record/5762297/files/student_download.zip?download=1)

## 👩‍🏫👩‍💻🎓 Instructors
Johan Hidding, Hanno Spreeuw

## 🧑‍🙋 Helpers
Dafne van Kuppevelt, Suvayu Ali 

## 👩‍💻👩‍💼👨‍🔬🧑‍🔬🧑‍🚀🧙‍♂️🔧 Roll Call
Name/ pronouns (optional) / job, role / social media (twitter, github, ...) / background or interests (optional) / city


## Icebreaker
Where did you experience the most beautiful night sky?


## 🗓️ Agenda

## 🔧 Exercises
### Exercise: getting meta-data
One of the other tables we will use is `gaiadr2.panstarrs1_original_valid`. Use `load_table` to get the metadata for this table. How many columns are there and what are their names?

Solution:
```python=
panstarrs_metadata = Gaia.load_table('gaiadr2.panstarrs1_original_valid')
print(panstarrs_metadata)

for column in panstarrs_metadata.columns:
    print(column.name)
```

### Exercise: read the * manual
Read the [documentation](https://gea.esac.esa.int/archive/documentation/GDR2/Gaia_archive/chap_datamodel/sec_dm_main_tables/ssec_dm_gaia_source.html) of this table and choose a column that looks interesting to you. Add the column name to the query and run it again. What are the units of the column you selected? What is its data type?


### Exercise: fix the query
The clauses in a query have to be in the right order. Go back and change the order of the clauses in `query2` and run it again. The modified query should fail, but notice that you don’t get much useful debugging information.

For this reason, developing and debugging ADQL queries can be really hard. A few suggestions that might help:

1. Whenever possible, start with a working query, either an example you find online or a query you have used in the past.
2. Make small changes and test each change before you continue.
3. While you are debugging, use TOP to limit the number of rows in the result. That will make each test run faster, which reduces your development time.
4. Launching test queries synchronously might make them start faster, too.


Solution:
Unfortunately, whatever you do wrong, you always get the same error! There is no way to get more detailed information.
Note that `SELECT ... FROM ... WHERE ...` always need to be in the same order.

### Exercise: operators
Read about SQL operators [here](https://www.w3schools.com/sql/sql_operators.asp) and then modify the previous query to select rows where `bp_rp` is between -0.75 and 2.


Solution:
Two options:
```
"""SELECT 
TOP 10
source_id, ra, dec, parallax, bp_rp
FROM gaiadr2.gaia_source
WHERE parallax < 1 
  AND bp_rp BETWEEN -0.75 AND 2
"""
```
Or:
```
"""SELECT 
TOP 10
source_id, ra, dec, parallax, bp_rp
FROM gaiadr2.gaia_source
WHERE parallax < 1 
  AND bp_rp > -0.75 AND bp_rp < 2
"""
```

### Exercise:
This query always selects sources with parallax less than 1. But suppose you want to take that upper bound as an input.

Modify `query3_base` to replace 1 with a format specifier like `{max_parallax}`. Now, when you call format, add a keyword argument that assigns a value to `max_parallax`, and confirm that the format specifier gets replaced with the value you provide.


```

Solution:
```python=
query3_base = """SELECT
TOP 10
{columns}
FROM gaiadr2.gaia_source
WHERE parallax < {max_parallax}
  AND bp_rp BETWEEN -0.75 AND 2
"""

query3 = query3_base.format(columns=','.join(columns), max_parallax=3)
```

### Exercise: arcminutes
Create a quantity that represents 5 arcminutes and assign it to a variable called radius.

Then convert it to degrees

Solution:
```python=
small_angle = 5 * u.arcmin
small_angle.to(u.deg)
```

### Exercise: TOP
When you are debugging queries like this, you can use `TOP` to limit the size of the results, but then you still don’t know how big the results will be.

An alternative is to use `COUNT`, which asks for the number of rows that would be selected, but it does not return them.

In the previous query, replace `TOP 10` source_id with `COUNT(source_id)` and run the query again. How many stars has Gaia identified in the cone we searched?

Solution:
```python=
count_cone_query = """SELECT 
COUNT(source_id)
FROM gaiadr2.gaia_source
WHERE 1=CONTAINS(
  POINT(ra, dec),
  CIRCLE(88.8, 7.4, 0.08333333))
"""

count_cone_job = Gaia.launch_job(count_cone_query)
count_cone_results = count_cone_job.get_results()
count_cone_results
```
(Note that there are other aggregate functions than `COUNT` available in ADQL)

### Exercise: ICRS
Find the location of GD-1 in ICRS coordinates.

- Create a SkyCoord object at 0°, 0° in the GD-1 frame.
- Transform it to the ICRS frame.

Hint: Because ICRS is a standard frame, it is built into Astropy. You can specify it by name, icrs (as we did with galactic).


```python=
gd1coord = SkyCoord(phi1=0*u.deg, phi2=0*u.deg, frame=gd1_frame)
gd1coord.transform_to('icrs')
```
<SkyCoord (ICRS): (ra, dec) in deg
    (200., 59.4504341)>
    
    
Solution:
```python=
SkyCoord(0*u.deg, 0*u.deg, frame=gd1_frame).transform('icrs')
```
Note that gd1 frame is defined with phi1 and phi2, not ra and dec.
Note also that you *need* to include units, otherwise you get an error!

### Exercise plotting
In the call to plt.plot, use the keyword argument `markersize` to make the markers smaller.

Then add the keyword argument `alpha` to make the markers partly transparent.

Adjust these arguments until you think the figure shows the data most clearly.


Solution:
```python=
plt.plot(x, y, 'ko', ms=0.2, alpha=0.1)

plt.xlabel('ra (degree ICRS)')
plt.ylabel('dec (degree ICRS)')
```

## 🧠 Collaborative Notes
If you followed the instructions to install conda, you can type in your terminal:
```
conda activate AstronomicalData
```
And make sure to install jupyter-lab:
```
conda install jupyterlab
```

and start:
```
jupyter lab
```

Make sure to navigate to a folder with student_data and create a new notebook.

Before we start:
* Do you have experience with Jupyter?
* Do you have jupyter ready to go?

Answers:

[shortkeys for Jupyter](https://towardsdatascience.com/jypyter-notebook-shortcuts-bf0101a98330)


```python=
# importing Gaia:
from astroquery.gaia import Gaia

# Connect to server to load names of tables:
tables = Gaia.load_tables(only_names=True)

# Print names with loop
for table in tables:
    print(table.name)
    
```
We are going to look at:
- `gaiadr2.gaia_source`
- `gaiadr2panstarrs1_original_valid`
- `gaiadr2.panstarrs1_best_neighbour`

```python=
# Load metadata of one table
table_metadata = Gaia.load_table('gaiadr2.gaia_source')
print(table.metadata)

for column in table_metadata.columns:
    print(column.name)
```

Let's write a query. We only take the first 10 rows, because we want to see if the query is working. It could be that we have millions of rows!
Queries are always of the form `SELECT ... FROM ...`
```python=
query1 = """SELECT
TOP 10
source_id, ra, parallax
FROM gaiadr2.gaia_source
"""

# Send to server:
job1 = Gaia.launch_job(query1)
print(job1)

# Retrieve results
results1 = job.get_results()
type(results1) # This is an astropy table
results1
```

Astropy tables are a bit similar to pandas tables, and are nicely shown in jupyter.

If we move to larger datasets, we need to use *asynchonous queries*. Synchronous queries are limited to 2000 rows.

We add a third keyword to our query: `SELECT ... FROM ... WHERE ...`
```python=
query2 = """SELECT
source_id, ra, dec, pmra, pmdec, parallax
TOP 3000
FROM gaiadr2.gaia_source
WHERE parallax < 1
"""

job2 = Gaia.launch_job_async(query2)
job2

results2 = job2.get_results()
results2
```
Results could be different from Johans, because order of records is not fixed. 

What is bp_rp? see this image:
![](https://datacarpentry.org/astronomy-python/fig/Gaia-HR-diagram.jpeg)


Formatting queries:
```python=
columns = 'source_id, ra, dec, pmra, pmdec, parallax'

query3_base = """SELECT
TOP 10
{columns}
FROM gaiadr2.gaia_source
WHERE parallax < 1
  AND bp_rp BETWEEN -0.75 AND 2
"""

query3 = query3_base.format(columns=columns)
```

Or what also works:
```python=
columns = ['source_id', 'ra', 'dec', 'pmra', 'pmdec', 'parallax']
query3_base = """SELECT 
TOP 10
{columns}
FROM gaiadr2.gaia_source
WHERE parallax < 1
  AND bp_rp BETWEEN -0.75 AND 2
"""

query3 = query3_base.format(columns=','.join(columns))
```

Let's look at the results:
```python=
job3 = Gaia.launch_job(query3)
print(job3)
results3 = job3.get_results()
results3

```

### Coordinate tranformations
(you can create a new notebook for this part)

```python=
import astropy.units as u

# Show all functionalities
dir(u)

angle = 10 * u.degree
type(angle)

# Convert
angle.to(u.arcmin)

# Can add quantities with matchable units
angle + 30*u.arcmin
```

We are going to select a cone region from the sky
```python=
cone_query = """SELECT
TOP 10
source_id
FROM gaiadr2.gaia_source
WHERE 1 = CONTAINS(
    POINT(ra,dec),
    CIRCLE(88.8, 7.4, 0.08333333))
"""
```
Note the new keywords `CONTAIN`, `POINT` and `CIRCLE` (we could have used other regions here such as polygons). They are specific to ADQL. The [CONTAINS function](https://www.ivoa.net/documents/ADQL/20180112/PR-ADQL-2.1-20180112.html#tth_sEc4.2.12) returns either 0 (when False) or 1 (when True). 
The points are all in ICRS coordinates, will need to do some conversion.

```python=
from astroquery.gaia import Gaia

cone_job = Gaia.launch_job(cone_query)
cone_job

cone_job.get_results()
```

Now we want to try to reproduce this image:
![](https://datacarpentry.org/astronomy-python/fig/gd1-4.png)
They work with their own coordinate system, so we need to convert.

We create a coordinates object `coord_icrs`.

```python=
from astropy.coordinates import SkyCoord

ra = 88.8 * u.degree
dec = 7.4 * u.degree

coord_icrs = SkyCoord(ra=ra, dec=dec, frame='icrs')

coord_icrs

# Can transform to different coordinate system:
coord_icrs.transform_to('galactic')
```

We need a separate library that contains a lot of coordinate systems to get the one from the paper
```python=
from gala.coordinates import GD1Koposov10

gd1_frame = GD1Koposov10()
type(gd1_frame)

# And transform to this frame:
coord_icrs.transform_to(gd1_frame)
```

Now we want to select a rectangle in the sky. We can not use the [BOX selector](https://www.ivoa.net/documents/ADQL/20180112/PR-ADQL-2.1-20180112.html#tth_sEc4.2.9) because we can only select a box in ICRS but we want in GD1! So we need to take the corners of the box and convert it to a polygon.

These are the corners in GD1:
```python=
phi1_min = -55 * u.degree
phi1_max = -45 * u.degree
phi2_min = -8 * u.degree
phi2_max = 4 * u.degree
```

We make a function for creating a rectangle with these corners (to prevent repetition):
```python=
def make_rectangle(x1, x2, y1, y2):
    """Return the corners of a rectangle."""
    xs = [x1, x1, x2, x2, x1] # coordinates in phi1
    ys = [y1, y2, y2, y1, y1] # coordinates in phi2
    return xs, ys

# And use it on our box:
phi1_rect, phi2_rect = make_rectangle(
    phi1_min, phi1_max, phi2_min, phi2_max)

# Now we can make a skycoord object with the corners in GD1:
corners = SkyCoord(phi1=phi1_rect, phi2=phi2_rect, frame=gd1_frame)
corners

# And tranform to ICRS:
corners_icrs = corners.transform_to('icrs')
corners_icrs
```

Now we can use the [`POLYGON` function](https://www.ivoa.net/documents/ADQL/20180112/PR-ADQL-2.1-20180112.html#tth_sEc4.2.19) 

To get our corners into the function, first need to convert to string:
```python=
def skycoord_to_string(skycord):
    """Convert a 1-d list of SkyCoord to string in ADQL format"""
    return ' '.join(skycord.to_string()).replace(' ', ',')

sky_point_list = skycoord_to_string(corners_icrs)


# put it inside the query:
columns = ['source_id', 'ra', 'dec', 'pmra', 'pmdec', 'parallax']
query3_base = """SELECT 
TOP 10 
{columns}
FROM gaiadr2.gaia_source
WHERE parallax < {max_parallax}
  AND bp_rp BETWEEN -0.75 AND 2 
  AND 1 = CONTAINS(POINT(ra, dec), 
                   POLYGON({sky_point_list}))
"""

query3 = query3_base.format(
    columns=','.join(columns),
    max_parallax=1,
    sky_point_list=sky_point_list)

job = Gaia.launch_job_async(query3)
job.get_results()
```

Now let's get the comple result (not top 10) and save the result:
```python=
columns = ['source_id', 'ra', 'dec', 'pmra', 'pmdec', 'parallax']
query3_base = """SELECT 
{columns}
FROM gaiadr2.gaia_source
WHERE parallax < {max_parallax}
  AND bp_rp BETWEEN -0.75 AND 2 
  AND 1 = CONTAINS(POINT(ra, dec), 
                   POLYGON({sky_point_list}))
"""

query3 = query3_base.format(
    columns=','.join(columns),
    max_parallax=1,
    sky_point_list=sky_point_list)
print(query3) # To check

job = Gaia.launch_job_async(query3)

# Save the results
polygon_results = job.get_results()
len(polygon_results)

filename = 'gd1_results.fits'
polygon_results.write(filename, overwrite=True)
```
This writes a file of 6,5MB. (there's a backup of this file in the `student_download` directory)

The complete query should look like this:

```sql=
SELECT
source_id,ra,dec,pmra,pmdec,parallax
FROM gaiadr2.gaia_source
WHERE parallax < 1
   AND bp_rp BETWEEN -0.75 AND 2
   AND 1 = CONTAINS(POINT(ra, dec),
               POLYGON(146.275,19.2619,135.422,25.8774,141.603,34.3048,152.817,27.1361))
```

### Plotting
(start a new notebook, need to do this to catch up:)
```python
from astroquery.gaia import Gaia
from astropy.table import Table

# read file from previous part
filename = 'gd1_results.fits' # 'backup-data/gd1_results.fits'
polygon_results = Table.read(filename)

# Create coordinate frame
from gala.coordinates import GD1Koposov10
gd1_frame = GD1Koposov10()
```

Look into the table that we read in:
```python=
polygon_results.info()

polygon_results.colnames

# Retrieve first row:
polygon_results[0]

# Select one column
polygon_results['ra']

# Or combine to get one value:
polygon_results[0]['ra']
```

For plotting, we need to choose x and y axis.
```python=
x = polygon_results['ra']
y = polygon_results['dec']

# Import the plotting library
import matplotlib.pyplot as plt

# Plot the results
plt.plot(x, y, 'ko') # ko stands for black dots
plt.xlabel('ra (degree ICRS)')
plt.ylabel('dec (degree ICRS)')
```

If this doesn't show anything, add after importing matplotlib
```
%matplotlib inline
```

Now we want to plot in the GD-1 frame, so we need to transform back.
```python=
from astropy.coordinates import SkyCoord
import astropy.units as u

distance = 8 * u.kpc
radial_velocity = 0 * u.km / u.sec

# Select everything that we're interested in 
skycoord = SkyCoord(ra=polygon_results['ra'], 
                    dec=polygon_results['dec'],
                    pm_ra_cosdec = polygon_results['pmra'],
                    pm_dec = polygon_results['pmdec'],
                    distance = distance,
                    radial_velocity = radial_velocity
                   )
```

Tomorrow we will continue with Reflex correction!

## Questions and answers
### SSL issues when downloading data from Gaia

For anyone encountering this error: "SSLCertVerificationError: [SSL: CERTIFICATE_VERIFY_FAILED] certificate verify failed: self signed certificate in certificate chain (_ssl.c:1056) "
We found a solution through [this link](https://community.netapp.com/t5/Software-Development-Kit-SDK-and-API-Discussions/Python-How-to-disable-SSL-certificate-verification/td-p/113697):
```
import ssl

try:
    _create_unverified_https_context = ssl._create_unverified_context
except AttributeError:
    # Legacy Python that doesn't verify HTTPS certificates by default
    pass
else:
    # Handle target environment that doesn't support HTTPS verification
    ssl._create_default_https_context = _create_unverified_https_context
```




### Difference between SQL and ADQL?
There are many flavors/dialects of SQL, ADQL is one of them but with the addition to work with astronomical coordinate system.

### What is Gaia?
Gaia is a satelite, measures position (and velocity) of stars with high precision

### Is u.deg same as u.degree?
Yes, you can check:
```python=
u.deg == u.degree
```

### Are there no units in the query?
No, but the documentation will tell you what the units are

## Tips and Tops
### Tip
Johan and Hanno had drastically different levels of preparedness, it was a difficult switch.

The last program of Johan was not working for me ... error in collaboration sheet?

Because the code is also put on the document live, due to tiny errors in writing / slight delay in typing it is sometimes hard to keep up with what the instructor is writing - for example when switching quickly to an exercise where you need the code in the notebook

Maybe it would be nice to have examples on how the things we see can be generalized in order to use them also for other types of data/catalogues.

### Top
Good preparation of the course. Everything works in the initial setup.

Plotting is so satisfying!

Coding live with the instructors is a nice way of being actively involved with the material 

Very nice to do hands-on exercises right from the beginning.

## 📚 Resources
- [NL-RSE meetup](https://www.eventbrite.co.uk/e/rse-career-stages-from-junior-to-group-lead-tickets-262466312807)
- [eScience Center Fellowship Programme](https://www.esciencecenter.nl/fellowship-programme/)
- [About the Table Access Protocol of astro data](https://ivoa.net/documents/TAP/)
- [ADQL documentation](https://www.ivoa.net/documents/ADQL/20180112/PR-ADQL-2.1-20180112.html)
- [AstroPy Units](https://docs.astropy.org/en/stable/units/)
- [The paper we are reproducing](https://iopscience.iop.org/article/10.3847/2041-8213/aad7b5)