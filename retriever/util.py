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


def make_location(start, end=None, strand=None):
    if end is not None:
        location = {'type': 'range',
                    'start': {'type': 'point',
                              'position': int(start)},
                    'end': {'type': 'point',
                            'position': int(end)}}
    else:
        location = {'type': 'point',
                    'position': int(start)},
    if strand is not None:
        location['strand'] = int(strand)
    return location
