#!/usr/bin/env python3
"""Replace generated slides with hand-edited versions stored in overrides/.

build_deck.js regenerates component_entropy_results.pptx from scratch, which
would discard edits made in PowerPoint. Each overrides/slideN/ directory holds
the hand-edited slide XML (slide.xml), its relationships (slide.xml.rels) and
any media it references (eqK.png). This script copies them into slide N of the
generated deck. Run it after `node build_deck.js`.

overrides/slide4/ keeps the native PowerPoint equations of the HNC slide
(commit df65b4d "corrected equations"). The override must reference the same
chart, layout and notes parts as the generated slide; this is checked.
"""
import re
import shutil
import sys
import tempfile
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
DECK = HERE / 'component_entropy_results.pptx'


def norm(target):
    return target.replace('/ppt/', '../', 1) if target.startswith('/ppt/') else target


def main(deck=DECK):
    with zipfile.ZipFile(deck) as source:
        parts = {name: source.read(name) for name in source.namelist()}
    for folder in sorted((HERE / 'overrides').glob('slide*')):
        number = int(folder.name[5:])
        slide = f'ppt/slides/slide{number}.xml'
        rels = f'ppt/slides/_rels/slide{number}.xml.rels'
        assert slide in parts, f'{slide} missing from generated deck'
        old_rels = parts[rels].decode()
        new_rels = (folder / 'slide.xml.rels').read_text()
        kinds = ('chart', 'slideLayout', 'notesSlide')
        for kind in kinds:
            old = {norm(t) for t in re.findall(rf'/{kind}" Target="([^"]+)"', old_rels)}
            new = {norm(t) for t in re.findall(rf'/{kind}" Target="([^"]+)"', new_rels)}
            assert old == new, f'slide {number}: {kind} targets differ ({old} vs {new})'
        images = sorted(folder.glob('eq*.png'))
        targets = re.findall(r'/image" Target="\.\./media/([^"]+)"', new_rels)
        assert len(targets) == len(images), f'slide {number}: image count mismatch'
        for index, target in enumerate(targets):
            name = f'override-{number}-{index + 1}.png'
            new_rels = new_rels.replace(f'../media/{target}"', f'../media/{name}"')
            parts[f'ppt/media/{name}'] = (folder / f'eq{target_index(target)}.png').read_bytes()
        parts[slide] = (folder / 'slide.xml').read_bytes()
        parts[rels] = new_rels.encode()
    types = parts['[Content_Types].xml'].decode()
    assert 'Extension="png"' in types
    with tempfile.NamedTemporaryFile(dir=deck.parent, suffix='.pptx', delete=False) as handle:
        tmp = Path(handle.name)
    with zipfile.ZipFile(tmp, 'w', zipfile.ZIP_DEFLATED) as target:
        if '[Content_Types].xml' in parts:
            target.writestr('[Content_Types].xml', parts.pop('[Content_Types].xml'))
        for name, data in parts.items():
            target.writestr(name, data)
    shutil.move(tmp, deck)
    print(f'Applied overrides to {deck.name}')


def target_index(target):
    """Map the PowerPoint media name imageK.png to the stored eqK.png."""
    return int(re.findall(r'\d+', target)[0])


if __name__ == '__main__':
    main(Path(sys.argv[1]) if len(sys.argv) > 1 else DECK)
