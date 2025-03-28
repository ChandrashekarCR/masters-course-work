import streamlit as st
import pandas as pd
import folium
from streamlit_folium import folium_static
from chromosome_plot import plot_chromosome_segments, plot_chromosome_country, plot_chromosome_continent
from world_map import plot_world_map, plot_world_map_country, plot_world_map_continent


# Hardcoded test data - Some random values. Use the test data for better results.

test_data = {'SAMPLE_ID': {0: 'HGDP00660',
  1: 'HGDP00660',
  2: 'HGDP00660',
  3: 'HGDP00660',
  4: 'HGDP00660',
  5: 'HGDP00660',
  6: 'HGDP00660',
  7: 'HGDP00660',
  8: 'HGDP00660',
  9: 'HGDP00660',
  10: 'HGDP00660',
  11: 'HGDP00660',
  12: 'HGDP00660',
  13: 'HGDP00660',
  14: 'HGDP00660',
  15: 'HGDP00660',
  16: 'HGDP00660',
  17: 'HGDP00660',
  18: 'HGDP00660',
  19: 'HGDP00660',
  20: 'HGDP00660',
  21: 'HGDP00660',
  22: 'HGDP00660',
  23: 'HGDP00660',
  24: 'HGDP00660',
  25: 'HGDP00660',
  26: 'HGDP00660'},
 'individual': {0: 351,
  1: 351,
  2: 351,
  3: 351,
  4: 351,
  5: 351,
  6: 351,
  7: 351,
  8: 351,
  9: 351,
  10: 351,
  11: 351,
  12: 351,
  13: 351,
  14: 351,
  15: 351,
  16: 351,
  17: 351,
  18: 351,
  19: 351,
  20: 351,
  21: 351,
  22: 351,
  23: 351,
  24: 351,
  25: 351,
  26: 351},
 'chromosome': {0: 1,
  1: 1,
  2: 1,
  3: 1,
  4: 1,
  5: 1,
  6: 2,
  7: 2,
  8: 2,
  9: 3,
  10: 3,
  11: 3,
  12: 3,
  13: 3,
  14: 3,
  15: 4,
  16: 4,
  17: 4,
  18: 4,
  19: 5,
  20: 5,
  21: 6,
  22: 7,
  23: 8,
  24: 9,
  25: 9,
  26: 9},
 'start_pos': {0: 1061166,
  1: 22638050,
  2: 50757885,
  3: 100051668,
  4: 185410975,
  5: 211785065,
  6: 12589208,
  7: 103176411,
  8: 131063169,
  9: 3846526,
  10: 22007867,
  11: 46310921,
  12: 71941154,
  13: 112072833,
  14: 169256240,
  15: 18473310,
  16: 44633416,
  17: 122744067,
  18: 156173215,
  19: 15571736,
  20: 152960378,
  21: 14329647,
  22: 3609712,
  23: 9375944,
  24: 1355061,
  25: 17626616,
  26: 75206337},
 'end_pos': {0: 22629057,
  1: 50722989,
  2: 100050789,
  3: 185344764,
  4: 211776439,
  5: 245770642,
  6: 103123301,
  7: 130665622,
  8: 229406167,
  9: 21996586,
  10: 46150400,
  11: 71930936,
  12: 112066141,
  13: 169212367,
  14: 193722330,
  15: 44531665,
  16: 122659704,
  17: 156091935,
  18: 185788762,
  19: 152940728,
  20: 175059713,
  21: 146960353,
  22: 146984783,
  23: 126554039,
  24: 17610211,
  25: 75201354,
  26: 120088933},
 'Prediction': {0: 'Tatars_2',
  1: 'Altaians_4',
  2: 'QUECHUA_3',
  3: 'Hidalgo_Mexico_1',
  4: 'Tatars_2',
  5: 'Hidalgo_Mexico_1',
  6: 'Hidalgo_Mexico_1',
  7: 'Mongols_2',
  8: 'Altaians_4',
  9: 'Altaians_3',
  10: 'Mongols_4',
  11: 'Altaians_4',
  12: 'Altaians_4',
  13: 'Mongols_4',
  14: 'Mongols_4',
  15: 'Altaians_3',
  16: 'Altaians_3',
  17: 'Mongols_4',
  18: 'Mongols_2',
  19: 'QUECHUA_1',
  20: 'Mongols_3',
  21: 'Hidalgo_Mexico_1',
  22: 'Hidalgo_Mexico_1',
  23: 'Hidalgo_Mexico_1',
  24: 'Altaians_3',
  25: 'Bermudian_2',
  26: 'Altaians_4'},
 'Population': {0: 'Surui',
  1: 'Surui',
  2: 'Surui',
  3: 'Surui',
  4: 'Surui',
  5: 'Surui',
  6: 'Surui',
  7: 'Surui',
  8: 'Surui',
  9: 'Surui',
  10: 'Surui',
  11: 'Surui',
  12: 'Surui',
  13: 'Surui',
  14: 'Surui',
  15: 'Surui',
  16: 'Surui',
  17: 'Surui',
  18: 'Surui',
  19: 'Surui',
  20: 'Surui',
  21: 'Surui',
  22: 'Surui',
  23: 'Surui',
  24: 'Surui',
  25: 'Surui',
  26: 'Surui'},
 'Lat': {0: 56.2901712697932,
  1: 49.5343807782193,
  2: -13.996779447178,
  3: 16.2724613000608,
  4: 54.3087225856653,
  5: 14.5388286401909,
  6: 14.5388286401909,
  7: 43.15549508082,
  8: 45.6742246652587,
  9: 48.1165424910612,
  10: 45.7297961757076,
  11: 43.8306322133876,
  12: 51.8994887304555,
  13: 45.340660487638,
  14: 44.7527713029204,
  15: 47.8650246078907,
  16: 46.6121418977288,
  17: 44.6572323476927,
  18: 43.0773796803656,
  19: -11.7584965128756,
  20: 43.3360833877639,
  21: 16.1283181828405,
  22: 10.134885565629,
  23: 14.1262181665565,
  24: 47.021514517507,
  25: 5.40985097902158,
  26: 48.7047229654416},
 'Lon': {0: 51.5119901277749,
  1: 86.7622967703801,
  2: -73.4656378403744,
  3: -95.7050004279268,
  4: 29.8005629616605,
  5: -92.2277500068698,
  6: -92.2277500068698,
  7: 115.579420750641,
  8: 90.9976028763527,
  9: 96.8374097054337,
  10: 106.352468649561,
  11: 88.4819737439493,
  12: 65.5041590679354,
  13: 107.393958004197,
  14: 110.025979685534,
  15: 93.9773206400829,
  16: 100.347063189851,
  17: 110.611552360222,
  18: 115.390955808296,
  19: -75.076996571367,
  20: 115.297918024514,
  21: -95.2502270169731,
  22: -85.7974448310629,
  23: -91.6897466702791,
  24: 100.951844691863,
  25: -52.8821412827541,
  26: 90.6323113192822},
 'Country': {0: 'Russian Federation',
  1: 'Russian Federation',
  2: 'Peru',
  3: 'Mexico',
  4: 'Belarus',
  5: 'Mexico',
  6: 'Mexico',
  7: 'China',
  8: 'China',
  9: 'Mongolia',
  10: 'Mongolia',
  11: 'China',
  12: 'Kazakhstan',
  13: 'Mongolia',
  14: 'Mongolia',
  15: 'Mongolia',
  16: 'Mongolia',
  17: 'Mongolia',
  18: 'China',
  19: 'Peru',
  20: 'China',
  21: 'Mexico',
  22: 'Costa Rica',
  23: 'Guatemala',
  24: 'Mongolia',
  25: 'French Guiana',
  26: 'Mongolia'},
 'Continent': {0: 'Europe',
  1: 'Europe',
  2: 'South America',
  3: 'North America',
  4: 'Europe',
  5: 'North America',
  6: 'North America',
  7: 'Asia',
  8: 'Asia',
  9: 'Asia',
  10: 'Asia',
  11: 'Asia',
  12: 'Asia',
  13: 'Asia',
  14: 'Asia',
  15: 'Asia',
  16: 'Asia',
  17: 'Asia',
  18: 'Asia',
  19: 'South America',
  20: 'Asia',
  21: 'North America',
  22: 'North America',
  23: 'North America',
  24: 'Asia',
  25: 'South America',
  26: 'Asia'},
 'Admixture1': {0: 1e-05,
  1: 1e-05,
  2: 1e-05,
  3: 1e-05,
  4: 1e-05,
  5: 1e-05,
  6: 0.00088475,
  7: 1e-05,
  8: 1e-05,
  9: 0.016057,
  10: 1e-05,
  11: 1e-05,
  12: 1e-05,
  13: 0.069389,
  14: 1e-05,
  15: 1e-05,
  16: 1e-05,
  17: 1e-05,
  18: 1e-05,
  19: 1e-05,
  20: 1e-05,
  21: 0.0149915,
  22: 1e-05,
  23: 1e-05,
  24: 1e-05,
  25: 1e-05,
  26: 0.0216665},
 'Admixture2': {0: 1e-05,
  1: 1e-05,
  2: 1e-05,
  3: 0.0315385,
  4: 1e-05,
  5: 1e-05,
  6: 1e-05,
  7: 0.092283,
  8: 0.015574,
  9: 1e-05,
  10: 0.021809,
  11: 1e-05,
  12: 1e-05,
  13: 0.0341815,
  14: 0.016578,
  15: 1e-05,
  16: 1e-05,
  17: 1e-05,
  18: 1e-05,
  19: 0.0444735,
  20: 1e-05,
  21: 0.0256072499999999,
  22: 0.00099175,
  23: 0.01747975,
  24: 1e-05,
  25: 1e-05,
  26: 0.0358099999999999},
 'Admixture3': {0: 1e-05,
  1: 0.096565,
  2: 1e-05,
  3: 1e-05,
  4: 1e-05,
  5: 0.062113,
  6: 0.0642119999999999,
  7: 0.012224,
  8: 0.05013375,
  9: 1e-05,
  10: 1e-05,
  11: 1e-05,
  12: 0.136359,
  13: 1e-05,
  14: 1e-05,
  15: 1e-05,
  16: 1e-05,
  17: 1e-05,
  18: 1e-05,
  19: 0.019697,
  20: 1e-05,
  21: 0.035462,
  22: 0.0890937499999999,
  23: 0.0040234999999999,
  24: 1e-05,
  25: 1e-05,
  26: 0.057762},
 'Admixture4': {0: 1e-05,
  1: 1e-05,
  2: 1e-05,
  3: 0.005795,
  4: 1e-05,
  5: 1e-05,
  6: 0.0178955,
  7: 1e-05,
  8: 0.0193072499999999,
  9: 1e-05,
  10: 1e-05,
  11: 1e-05,
  12: 1e-05,
  13: 0.0001025,
  14: 1e-05,
  15: 1e-05,
  16: 0.067177,
  17: 1e-05,
  18: 1.5e-05,
  19: 1.5e-05,
  20: 0.064682,
  21: 0.0016354999999999,
  22: 0.010071,
  23: 1e-05,
  24: 1e-05,
  25: 1e-05,
  26: 1e-05},
 'Admixture5': {0: 0.925113,
  1: 0.618949,
  2: 0.95179,
  3: 0.82096,
  4: 0.784826,
  5: 0.7430975,
  6: 0.80588075,
  7: 0.779821,
  8: 0.82498275,
  9: 0.772706,
  10: 0.808521,
  11: 0.89531,
  12: 0.662958,
  13: 0.771503,
  14: 0.906506,
  15: 0.792411,
  16: 0.913999,
  17: 0.918378,
  18: 0.868084,
  19: 0.856906,
  20: 0.822923,
  21: 0.805787,
  22: 0.69611325,
  23: 0.80803575,
  24: 0.939033,
  25: 0.888253,
  26: 0.7317415},
 'Admixture6': {0: 1e-05,
  1: 1e-05,
  2: 0.014834,
  3: 1e-05,
  4: 1e-05,
  5: 0.008943,
  6: 1e-05,
  7: 1e-05,
  8: 0.004482,
  9: 1e-05,
  10: 1e-05,
  11: 1e-05,
  12: 0.026124,
  13: 1e-05,
  14: 0.066335,
  15: 1e-05,
  16: 1e-05,
  17: 1e-05,
  18: 1e-05,
  19: 1e-05,
  20: 0.063668,
  21: 0.011314,
  22: 0.00355875,
  23: 0.0099225,
  24: 1e-05,
  25: 1e-05,
  26: 0.018945},
 'Admixture7': {0: 0.0499,
  1: 0.284426,
  2: 0.033316,
  3: 0.0755015,
  4: 0.215104,
  5: 0.107306,
  6: 0.11108675,
  7: 0.115623,
  8: 0.07779975,
  9: 0.211177,
  10: 0.019196,
  11: 0.075961,
  12: 0.129418,
  13: 0.1148515,
  14: 0.004474,
  15: 0.207519,
  16: 1e-05,
  17: 0.081552,
  18: 0.131841,
  19: 0.0661555,
  20: 0.048676,
  21: 0.08454225,
  22: 0.18319675,
  23: 0.156265,
  24: 0.060897,
  25: 0.111677,
  26: 0.134045},
 'Admixture8': {0: 0.024926,
  1: 1e-05,
  2: 1e-05,
  3: 0.0559105,
  4: 1e-05,
  5: 0.0785,
  6: 1e-05,
  7: 1e-05,
  8: 0.0023927499999999,
  9: 1e-05,
  10: 1e-05,
  11: 0.028669,
  12: 0.045101,
  13: 1e-05,
  14: 1e-05,
  15: 1e-05,
  16: 1e-05,
  17: 1e-05,
  18: 1e-05,
  19: 1e-05,
  20: 1e-05,
  21: 0.01379275,
  22: 0.0169549999999999,
  23: 1e-05,
  24: 1e-05,
  25: 1e-05,
  26: 1e-05},
 'Admixture9': {0: 1e-05,
  1: 1e-05,
  2: 1e-05,
  3: 0.010264,
  4: 1e-05,
  5: 1e-05,
  6: 1e-05,
  7: 1e-05,
  8: 0.0053177499999999,
  9: 1e-05,
  10: 0.150425,
  11: 1e-05,
  12: 1e-05,
  13: 0.009942,
  14: 0.006067,
  15: 1e-05,
  16: 0.018765,
  17: 1e-05,
  18: 1e-05,
  19: 0.0127225,
  20: 1e-05,
  21: 0.0068682499999999,
  22: 1e-05,
  23: 0.0042429999999999,
  24: 1e-05,
  25: 1e-05,
  26: 1e-05},
 'segment': {0: 1,
  1: 2,
  2: 3,
  3: 4,
  4: 5,
  5: 6,
  6: 1,
  7: 2,
  8: 3,
  9: 1,
  10: 2,
  11: 3,
  12: 4,
  13: 5,
  14: 6,
  15: 1,
  16: 2,
  17: 3,
  18: 4,
  19: 1,
  20: 2,
  21: 1,
  22: 1,
  23: 1,
  24: 1,
  25: 2,
  26: 3}}
country_mapping = {
    'Iran, Islamic Republic of': 'Iran',
    'Russian Federation': 'Russia',
    'Syrian Arab Republic': 'Syria',
    'Türkiye': 'Turkey',
    'Viet Nam': 'Vietnam',
    'Korea, Democratic People\'s Republic of': 'North Korea',
    'Korea, Republic of': 'South Korea',
    'Moldova, Republic of': 'Moldova',
    'Palestine, State of': 'Palestine',
    'Congo, The Democratic Republic of the': 'Dem. Rep. Congo',
    'Tanzania, United Republic of': 'Tanzania',
    'Lao People\'s Democratic Republic': 'Laos',
    'Venezuela, Bolivarian Republic of': 'Venezuela',
    'United States': 'United States of America',
    'Bolivia, Plurinational State of': 'Bolivia',
    'Brunei Darussalam': 'Brunei',
    'Cape Verde': 'Cabo Verde',
    'Czech Republic': 'Czechia',
    'Swaziland': 'eSwatini',
    'Gambia, The': 'Gambia',
    'Micronesia, Federated States of': 'Micronesia',
    'Burma': 'Myanmar',
    'The Bahamas': 'Bahamas',
    'The Gambia': 'Gambia',
    'Trinidad and Tobago': 'Trinidad and Tobago',
    'Falkland Islands (Malvinas)': 'Falkland Is.',
    'Timor-Leste': 'Timor-Leste',
    'Congo': 'Congo',
    'Ivory Coast': "Côte d'Ivoire",
    'Guinea-Bissau': 'Guinea-Bissau',
    'Eq. Guinea': 'Equatorial Guinea',
    'Central African Republic': 'Central African Rep.',
    'South Sudan': 'S. Sudan',
    'Western Sahara': 'W. Sahara',
    'Bosnia and Herzegovina': 'Bosnia and Herz.',
    'North Macedonia': 'North Macedonia',
    'United Arab Emirates': 'United Arab Emirates',
    'Saint Vincent and the Grenadines': 'Saint Vincent and the Grenadines',
    'Antigua and Barbuda': 'Antigua and Barbuda',
    'São Tomé and Príncipe': 'São Tomé and Principe',
    'Comoros': 'Comoros',
    'Saint Kitts and Nevis': 'Saint Kitts and Nevis',
    'Saint Lucia': 'Saint Lucia',
    'Seychelles': 'Seychelles',
    'Andorra': 'Andorra',
    'San Marino': 'San Marino',
    'Liechtenstein': 'Liechtenstein',
    'Monaco': 'Monaco',
    'Malta': 'Malta',
    'Marshall Islands': 'Marshall Islands',
    'Nauru': 'Nauru',
    'Palau': 'Palau',
    'Tuvalu': 'Tuvalu',
    'Vatican City': 'Vatican',
    'Solomon Islands': 'Solomon Is.',
    'New Zealand': 'New Zealand',
    'Guinea Bissau': 'Guinea-Bissau',
    'Dominican Republic': 'Dominican Rep.',
    'United Kingdom': 'United Kingdom',
    'Unknown': None  # Handle missing cases
}

# Convert test data to DataFrame
test_df = pd.DataFrame(test_data)

# Streamlit App Title
st.title("Ancestry Lens")

st.markdown("""
## About Ancestry Lens
Ancestry Lens is a visualization tool designed to explore **local ancestry** at the chromosomal level.  
Unlike global ancestry, which provides broad population estimates (e.g., 50% Italian, 50% English), local ancestry allows users to analyze **ancestry composition along different segments of chromosomes**.

### How to Use:
1. **Upload your CSV file** containing local ancestry data (must follow the correct format).  
   - If you don’t have data, a **test dataset** is available for exploration.  
   - This dataset includes a default population (**Surui**) and a pre-selected individual (**351**) with interactive chromosome selection.  
2. **Select a population** to filter individuals belonging to that group.  
3. **Select an individual** to visualize their chromosome segments.  
4. **Explore ancestry through various plots**:  
   - Chromosome-based visualizations (prediction, country, and continent levels).  
   - Interactive world maps highlighting ancestry distribution.  
5. **Refine the visualization by selecting a chromosome**, displaying its geographical distribution by country and continent.  
6. **Download all maps and plots as PNG files** for further analysis and sharing.  
7. **Ensure your file is in the correct format**:  
   - A **CSV file** with properly structured columns.  
   - You can check the example dataset at the bottom of the page and download a sample file for reference.  
8. **For users with Q files, BIM, FAM, and HGDP files**:  
   - Use the **preprocessing script** to convert these files into the correct CSV format before uploading.  

Enjoy discovering genetic histories with **Ancestry Lens**!  
""")


# Upload CSV File
uploaded_file = st.file_uploader("Upload CSV file", type=["csv"])

# Use test data if no file is uploaded
if uploaded_file:
    df = pd.read_csv(uploaded_file)
    df['Country'] = df['Country'].apply(lambda x: country_mapping.get(x,x))

else:
    df = test_df

# Sidebar filters
st.sidebar.header("Select Population and Individual")

# Select Population
selected_population = st.sidebar.selectbox("Select Population", ["None"] + list(df['Population'].dropna().unique()))

if selected_population != "None":
    filtered_df = df[df['Population'] == selected_population]
    individuals = filtered_df['individual'].unique()
else:
    filtered_df = df
    individuals = filtered_df['individual'].unique()

# Select Individual
selected_individual = st.sidebar.selectbox("Select Individual", ["None"] + list(individuals))

if selected_individual != "None":
    filtered_df = filtered_df[filtered_df['individual'] == selected_individual]
    filtered_df['segment'] = filtered_df.groupby('chromosome').cumcount() + 1

    # Chromosome Visualization
    with st.expander("Chromosome Visualization by Prediction"):
        st.plotly_chart(plot_chromosome_segments(filtered_df, selected_individual), use_container_width=True)

    with st.expander("Chromosome Visualization by Country"):
        st.plotly_chart(plot_chromosome_country(filtered_df, selected_individual), use_container_width=True)

    with st.expander("Chromosome Visualization by Continent"):
        st.plotly_chart(plot_chromosome_continent(filtered_df, selected_individual), use_container_width=True)

    # World Map
    with st.expander("World Map of Samples"):
        folium_static(plot_world_map(filtered_df, selected_individual))

    # Now ask user for chromosome to plot the world map by chromosome
    chromosomes = filtered_df['chromosome'].unique()
    selected_chromosome = st.sidebar.selectbox("Select Chromosome", chromosomes)

    if selected_chromosome:
        with st.expander(f"World Map of Chromosome {selected_chromosome} Distribution by Country"):
            st.plotly_chart(
                plot_world_map_country(
                    filtered_df, selected_chromosome, 
                    '/home/inf-21-2024/binp29/population_genetic_project/data/04_geopandas/ne_110m_admin_0_countries.shp'
                ), 
                use_container_width=True
            )

        with st.expander(f"World Map of Chromosome {selected_chromosome} Distribution by Continent"):
            st.plotly_chart(
                plot_world_map_continent(
                    filtered_df, selected_chromosome, 
                    '/home/inf-21-2024/binp29/population_genetic_project/data/04_geopandas/ne_110m_admin_0_countries.shp'
                ), 
                use_container_width=True
            )

    # Rename Admixture Columns # Make the change in country names here.
    admix_mapping = {
        'Admixture3': 'NORTHEASTASIAN', 'Admixture1': 'MEDITERRANIAN',
        'Admixture6': 'SOUTHAFRICA', 'Admixture8': 'SOUTHWESTASIAN',
        'Admixture2': 'NATIVEAMERICAN', 'Admixture5': 'OCEANIAN',
        'Admixture7': 'SOUTHEASTASIA', 'Admixture4': 'NORTHERNEUROPEAN',
        'Admixture9': 'SUBSAHARANAFRICA'
    }
    filtered_df.rename(columns=admix_mapping, inplace=True)

    # Expandable Full Data Table
    with st.expander(f"Full Data for Individual {selected_individual}"):
        st.dataframe(filtered_df, width=1000, height=500)