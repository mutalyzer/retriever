import requests


def make_request(url, params=None, headers=None):

    try:
        request = requests.get(url, params=params, headers=headers)
        request.raise_for_status()
    except requests.exceptions.HTTPError:
        return
    except requests.exceptions.ConnectionError as errc:
        print("Connection Error:", errc)
    except requests.exceptions.Timeout as errt:
        print("Timeout Error:", errt)
    except requests.exceptions.RequestException as err:
        print("Some other Error", err)
    else:
        return request.text
