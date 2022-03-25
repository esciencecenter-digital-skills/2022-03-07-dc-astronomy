![](https://i.imgur.com/iywjz8s.png)


# 2022-03-07 Astronomical Data Science (Day 3)

Welcome to The Workshop Collaborative Document.

This Document is synchronized as you type, so that everyone viewing this page sees the same text. This allows you to collaborate seamlessly on documents.

----------------------------------------------------------------------------

This is the Document for today: [link](https://hackmd.io/@O0tsDNPbTlyhyGiiCMaLIw/SykxkCS-9)

Collaborative Document day 1: [link](https://hackmd.io/@O0tsDNPbTlyhyGiiCMaLIw/SkeAMVmZ9)

Collaborative Document day 2: [link](https://hackmd.io/t25qWXHRRqWgvN6SVwcvAw)

Collaborative Document day 3: [link](https://hackmd.io/@O0tsDNPbTlyhyGiiCMaLIw/SykxkCS-9)

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
[link](<url>)

## 👩‍🏫👩‍💻🎓 Instructors
Hanno Spreeuw, Johan Hidding

## 🧑‍🙋 Helpers
Dafne van Kuppevelt, Suvayu Ali

## 👩‍💻👩‍💼👨‍🔬🧑‍🔬🧑‍🚀🧙‍♂️🔧 Roll Call
Name/ pronouns (optional) / job, role / social media (twitter, github, ...) / background or interests (optional) / city


## 🗓️ Agenda

| time | what |
|---|---|
| 9:00  |	Welcome and  Icebreaker |
| 9:15  |	JOIN  |
| 10:15 | 	break  |
| 10:30 |	More JOIN and Photometry |
| 11:30 | 	break  |
| 11:45 |	More Photometry and maybe a bit of Visualization |
| 12:45 | 	Closing and feedback  |

## 🔧 Exercises

### Check if you can run workaround


### Check if you have the HDF5 file from Teams


### Exercise: plot the 7345 stars
One dot per star along right ascension and declination
```python=
x = candidate_gaia_table["ra"]
y = candidate_gaia_table["dec"]

# plot commands
```


## 🧠 Collaborative Notes

### Following up on yesterday

There is a workaround for the `HTTPError` we got on our final query to the Gaia database:
- *reason:* proper motion values on Gaia are -ve (bug on their side)

**Note:** `pm_vertices` is same as day 2 [see collaborative notes](https://hackmd.io/t25qWXHRRqWgvN6SVwcvAw?view#%F0%9F%A7%A0-Collaborative-Notes).  If you are starting afresh, you could copy the following code to check if the workaround works for you:
```python=
import numpy as np
import pandas as pd

# you can also read from the HDF5 file you saved on day 2
selected_df = pd.read_hdf("student_download/data/"+filename, "selected_df")

from scipy.spatial import ConvexHull

hull  = ConvexHull(points)
pm_vertices = points[hull.vertices]

s = np.array2string(-pm_vertices.flatten(), max_line_width=1000, separator=",")
pm_point_list = s.strip("[]")

columns = "source_id, ra, dec, pmra, pmdec"
candidate_coord_query_base = """SELECT
{columns}
FROM gaiadr2.gaia_source
WHERE parallax < 1
  AND bp_rp BETWEEN -0.75 AND 2 
  AND 1 = CONTAINS(POINT(ra, dec), 
                   POLYGON({sky_point_list}))
  AND 1 = CONTAINS(POINT(-pmra, -pmdec),
                   POLYGON({pm_point_list}))
"""
# This a workaround for the bug on Gaia side, not recommended normally

candidate_coord_pm_query = candidate_coord_query_base.format(
    columns=columns,
    sky_point_list=sky_point_list,
    pm_point_list=pm_point_list
)

candidate_coord_pm_job = Gaia.launch_job_async(candidate_coord_pm_query)
print(candidate_coord_pm_job)

candidate_gaia_table = candidate_coord_pm_job.get_results()
len(candidate_gaia_table) # 7345
```

You could also read the results from [this HDF5 file](https://nlesc.sharepoint.com/:u:/s/instructors/EaWve28rtlZBiTRmoXb87uABRzgX7gLuLvNcCXSMeavr0g?e=wwRj1y):
```python
# note: this is a pandas.DataFrame, not an astropy table
candidate_gaia_table_df = pd.read_hdf("candidate_gaia_table.hdf", "candidate_gaia_table")
```

Let's do some plotting as a cross-check, [exercise](https://hackmd.io/oxqwu6EKRquC1pasy0W4qQ?both#Exercise-plot-the-7345-stars):
```python=
d = dict(sky_point_list=sky_point_list, pm_point_list=pm_point_list)
point_series = pd.Series(d)

# filename = "gd1_data.hdf"
point_series.to_hdf(filename, "point_series") # save checkpoint

x = candidate_gaia_table_df["ra"]
y = candidate_gaia_table_df["dec"]

plt.plot(x, y, 'ko', ms = 0.5, alpha=0.1)
plt.xlabel('ra(deg)')
plt.ylabel('dec(deg)')
plt.show()
```
We can import `make_dataframe` from before like this:
```python=
from episode_functions import *

type(make_dataframe)
```

If you don´t have `make_dataframe`, here it is:
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

### `JOIN` we will link different tables to get the full set of stars

```python=
ps_best_neighbor_meta = Gaia.load_table("gaiadr2.panstrarrs1_best_beighbour")
```

Gaia database is down, so we couldn't proceed further :frowning:

## Tips
(feedback for future iterations)

Perhaps send the notebooks after the day or before the next day. My pc went death and the notebooks were not entirely saved. The workaround did not work for me. The code needed more pre-definitions.

All the queries can be run in https://gaia.aip.de/query/ if the connection to Gaia is not working. If anyone like to follow the steps on the notes, they can do the queries on the website and download the results.

## Tops
(what you liked today)

Nice course. Let us continue on another day?

A great overview of how to query databases and get coordinates swapped around
Bad luck, but thank you for trying work-arounds!
I've checked the documentation step-by-step of the course, and that's very well done, gold! Great way to overcome problems like servers down...

## 📚 Resources

