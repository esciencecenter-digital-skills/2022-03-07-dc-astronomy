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
