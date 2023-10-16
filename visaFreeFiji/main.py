import requests
from bs4 import BeautifulSoup

# Specify the URL of the website to scrape
url = "https://visaguide.world/visa-free-countries/fijian-passport/"

# Send a GET request to the website and retrieve the HTML content
response = requests.get(url)
html_content = response.text

# Parse the HTML content using BeautifulSoup
soup = BeautifulSoup(html_content, "html.parser")

# Find the relevant section containing the visa-free countries
h2_header = soup.find("h2", {"id": "where-can-fijian-passport-holders-travel-without-a-visa"})
countries= h2_header.find


# Extract the list of countries
country_elements = countries_section.find_all("li")
visa_free_countries = [country.text.strip() for country in country_elements]

# Print the visa-free countries for Fiji passport holders
print("Visa-free countries for Fiji passport holders:")
for country in visa_free_countries:
    print(country)
