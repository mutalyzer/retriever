import requests


def get_gff(feature_id):
    url = 'https://rest.ensembl.org/lookup/id/{}'.format(feature_id)
    params = {'expand': 1,
              'Content-type': 'application/json'},

    try:
        response = requests.get(url, params=params)
        response.raise_for_status()
    except requests.exceptions.HTTPError as errh:
        print("HTTP Error:", errh)
    except requests.exceptions.ConnectionError as errc:
        print("Connection Error:", errc)
    except requests.exceptions.Timeout as errt:
        print("Timeout Error:", errt)
    except requests.exceptions.RequestException as err:
        print("Some other Error", err)
    else:
        print(response.headers.get('content-type'))
        return response.text
