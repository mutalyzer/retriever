import time

import requests


class RequestErrors(Exception):
    def __init__(self, errors):
        self.errors = errors


class Http400(Exception):
    def __init__(self, error):
        self.status_code = error.response.status_code
        self.response = error.response


def request(url, params=None, headers=None, timeout=1, max_retries=1, sleep=1):
    retries = 0
    errors = []
    while retries < max_retries:
        try:
            response = requests.get(
                url, params=params, headers=headers, timeout=timeout
            )
            response.raise_for_status()
        except requests.exceptions.HTTPError as e:
            if e.response.status_code == 400:
                raise Http400(e)
            errors.append(e)
            if e.response.status_code == 429:
                time.sleep(sleep)
        except Exception as e:
            errors.append(e)
        else:
            return response.text
        retries += 1

    raise RequestErrors(errors)
