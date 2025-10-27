# How to create the html for the github pages


Create the files `_config.yml` and `_toc.yml`. Then in the cmd realise:


```bash
python3 -m venv venv-book
```

```bash
source venv-book/bin/activate
```

```bash
source venv-book/bin/activate
```
```bash
pip install jupyter-book ghp-import
```

```bash
jupyter-book build .
```
In case of need, before rebuild: `jupyter-book clean .`

Then push to your github pages
```bash
ghp-import -n -p -f _build/html
```

Finally deactivate the environment:
```bash
deactivate
```

Add the `index.html` file to ensure starting in the good readme.