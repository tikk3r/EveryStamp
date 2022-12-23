import os
import sys
import logging

import tqdm
import requests

class FileDownloader(object):
    ''' From https://medium.com/better-programming/python-progress-bars-with-tqdm-by-example-ce98dbbc9697
    Copyright 2019 tiptapcode Authors. All Rights Reserved.
    Licensed under the Apache License, Version 2.0 (the "License");
    you may not use this file except in compliance with the License.
    You may obtain a copy of the License at
         http://www.apache.org/licenses/LICENSE-2.0
    Unless required by applicable law or agreed to in writing, software
    distributed under the License is distributed on an "AS IS" BASIS,
    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and
    limitations under the License.
    '''
    def get_url_filename(self, url):
        """
        Discover file name from HTTP URL, If none is discovered derive name from http redirect HTTP content header Location
        :param url: Url link to file to download
        :type url: str
        :return: Base filename
        :rtype: str
        """
        try:
            filename = os.path.basename(url)
            basename, ext = os.path.splitext(filename)
            if ext:
                return filename
            header = requests.head(url, allow_redirects=False).headers
            return os.path.basename(header.get('Location')) if 'Location' in header else filename
        except requests.exceptions.HTTPError as errh:
            print("Http Error:", errh)
            raise errh
        except requests.exceptions.ConnectionError as errc:
            print("Error Connecting:", errc)
            raise errc
        except requests.exceptions.Timeout as errt:
            print("Timeout Error:", errt)
            raise errt
        except requests.exceptions.RequestException as err:
            print("OOps: Something Else", err)
            raise err

    def download_file(self, url, filename=None, target_dir=None):
        """
        Stream downloads files via HTTP
        :param url: Url link to file to download
        :type url: str
        :param filename: filename overrides filename defined in Url param
        :type filename: str
        :param target_dir: target destination directory to download file to
        :type target_dir: str
        :return: Absolute path to target destination where file has been downloaded to
        :rtype: str
        """
        if target_dir and not os.path.isdir(target_dir):
            raise ValueError('Invalid target_dir={} specified'.format(target_dir))
        local_filename = self.get_url_filename(url) if not filename else filename

        req = requests.get(url, stream=True, verify=True)
        req.raise_for_status()
        file_size = int(req.headers['Content-Length'])
        chunk_size = 1024  # 1 MB
        num_bars = int(file_size / chunk_size)

        base_path = os.path.abspath(os.path.dirname(__file__))
        target_dest_dir = os.path.join(base_path, local_filename) if not target_dir else os.path.join(target_dir, local_filename)
        with open(target_dest_dir, 'wb') as fp:
            for chunk in tqdm.tqdm(req.iter_content(chunk_size=chunk_size), total=num_bars, unit='KB', leave=True, file=sys.stdout):
                fp.write(chunk)

        return target_dest_dir


class LegacyDownloader(FileDownloader):
    ''' Downloader sub-class for the DESI Legacy Imaging Surveys.
    '''
    supported_keywords = ['ra', 'dec', 'mode', 'layer', 'pixscale', 'bands', 'size_pix']
    logger = logging.getLogger('EveryStamp:LegacyDownloader')

    def __init__(self):
        self.url = 'https://www.legacysurvey.org/viewer/{mode:s}-cutout/?ra={ra:f}&dec={dec:f}&layer={layer:s}&pixscale={pixscale:.3f}&bands={bands:s}&size={size_pix:d}'

    def format_url(self, ra=None, dec=None, size=None, bands='grz', mode='jpeg', layer='ls-dr9', pixscale=0.262, autoscale=False, **kwargs):
        '''Returns a properly formatted URL that can be used to obtain a cutout from Legacy.
        '''
        size_pix = int(size * 3600 / pixscale)
        dlpixscale = pixscale
        dlsize_pix = size_pix
        if (size_pix > 3000) and autoscale:
            # Jump to the next available pixel size by scaling from the (approximate) native pixel scale.
            new_pixscale = pixscale
            new_size_pix = int(size * 3600 / new_pixscale)
            while new_size_pix > 3000:
                new_pixscale += 0.262
                new_size_pix = int(size * 3600 / new_pixscale)
            self.logger.warn('Image size of {:.2f} deg with pixel scale {:.3f} exceeds server limit of 3000 pixels! Automatically adjusting pixel scale to {:.3f} giving {:d} pixels.'.format(size, pixscale, new_pixscale, new_size_pix), Warning, stacklevel=2)
            dlpixscale = new_pixscale
            dlsize_pix = new_size_pix
        elif size_pix > 3000:
            self.logger.warn('Image size of {:.2f} deg with pixel scale {:.3f} exceeds server limit of 3000 pixels! Image will be truncated! Use --autoscale or pass autoscale=True to automatically switch pixel scales.'.format(size, pixscale), Warning, stacklevel=2)
        return self.url.format(ra=ra, dec=dec, size_pix=dlsize_pix, bands=bands, mode=mode, layer=layer, pixscale=dlpixscale)

    def download(self, **kwargs):
        furl = self.format_url(**kwargs)
        self.logger.info('Downloading cutout from %s', furl)
        if kwargs['ddir']:
            fname = kwargs['ddir'] + '/legacystamps_{ra:f}_{dec:f}_{layer:s}.{mode:s}'.format(ra=kwargs['ra'], dec=kwargs['dec'], layer=kwargs['layer'], mode=kwargs['mode'])
        else:
            self.logger.info('Download directory not specified, downloading to %s instead', os.getcwd())
        fname = os.getcwd() + '/legacystamps_{ra:f}_{dec:f}_{layer:s}.{mode:s}'.format(ra=kwargs['ra'], dec=kwargs['dec'], layer=kwargs['layer'], mode=kwargs['mode'])
        self.download_file(furl, filename=fname)
