---
layout: post
title:  "git-revert Vs. git-reset"
author: Rob
categories: [ Git ]
featured: true
---


# `$ git revert` Vs. `$ git reset`

### git-revert

- creates a new commit that undoes the changes from a previous commit. 
   - adds new history

```bash
± git add [file]; git commit -m "message"; git push

Enumerating objects:
...
   dc93c84..3dac56a  master -> master

# Note: 3dac56a is the parent of the latest commit.

± git revert 3dac56a 

[master 7321501] Revert "here goes the mistake..."
 1 file changed, 1 deletion(-)
 delete mode 100644 mistake.txt

± git push

# NOTE: you can check your commits with
± git rev-list --all --abbrev-commit HEAD
7321501
3dac56a
dc93c84

```




### git-reset

- Reset current HEAD to the specified state
   - could possibly remove history

```bash
± git add [file]; git commit -m "message"; git push

Enumerating objects:
...
   dc93c84..3dac56a  master -> master

# Note:  dc93c84 is the parent of the commit you wish to go back to.

± git reset --hard dc93c84

HEAD is now at dc93c84 Revert "message"

± git push -f




```






