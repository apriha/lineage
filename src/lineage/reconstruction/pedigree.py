"""The pedigree model: :class:`Person` and :class:`Pedigree`."""

from __future__ import annotations

import graphlib
import unicodedata
from collections.abc import Iterable, Iterator, Mapping
from dataclasses import dataclass, field, replace
from types import MappingProxyType
from typing import Self, TypeAlias

from lineage.exceptions import PedigreeError
from lineage.types import Sex

__all__ = ["Pedigree", "Person"]

_Record: TypeAlias = tuple[str, Sex, str | None, str | None]

_MAX_NAMES_SHOWN = 50


@dataclass(frozen=True, init=False, slots=True)
class Person:
    """A person in a pedigree.

    Parameters
    ----------
    name : str
        The person's name, unique within a pedigree. Names are case-sensitive and must
        not be empty, start or end with whitespace, or contain control characters such
        as tabs or newlines.
    sex : Sex or {"Male", "Female"}
        Stored as a :class:`~lineage.types.Sex`.
    mother, father : Person or str, optional
        The parent or the parent's name. Only the name is stored; a pedigree checks
        that a parent given as a ``Person`` matches the person it lists. Leave out a
        parent who is not in the pedigree.

    Raises
    ------
    TypeError
        If a name is not a str, or a parent is not a Person, a str or None.
    ValueError
        If a name is not valid, or `sex` is not "Male" or "Female".
    """

    name: str
    sex: Sex
    mother: str | None
    father: str | None
    _parents_given: tuple[_Record | None, _Record | None] = field(
        init=False, repr=False, compare=False
    )

    def __init__(
        self,
        name: str,
        sex: Sex | str,
        mother: Person | str | None = None,
        father: Person | str | None = None,
    ) -> None:
        name = _check_name(name, "A person's name")
        checked_sex = _check_sex(sex, name)
        mother_name, mother_given = _parent_link(mother, "mother", name)
        father_name, father_given = _parent_link(father, "father", name)
        object.__setattr__(self, "name", name)
        object.__setattr__(self, "sex", checked_sex)
        object.__setattr__(self, "mother", mother_name)
        object.__setattr__(self, "father", father_name)
        object.__setattr__(self, "_parents_given", (mother_given, father_given))


class Pedigree:
    """People and the parent links between them.

    A pedigree is immutable. Two pedigrees are equal when they hold the same people, in
    any order.

    Parameters
    ----------
    people : iterable of Person
        Everyone in the pedigree, in any order. Every parent a person names must be in
        it.

    Raises
    ------
    TypeError
        If an item of `people` is not a :class:`Person`.
    PedigreeError
        If two people have the same name; a person is their own parent, has the same
        person as mother and father, or has a parent who is not in the pedigree; a
        mother is not female or a father not male; a parent given as a
        :class:`Person` differs from the person listed under that name; or someone is
        their own ancestor.
    """

    __slots__ = (
        "_by_name",
        "_children",
        "_founders",
        "_hash",
        "_names",
        "_order",
        "_people",
    )

    def __init__(self, people: Iterable[Person]) -> None:
        items: tuple[object, ...] = tuple(people)
        by_name: dict[str, Person] = {}
        for item in items:
            if not isinstance(item, Person):
                raise TypeError(
                    f"Every item of a pedigree must be a Person, not "
                    f"{type(item).__name__}"
                )
            if item.name in by_name:
                raise PedigreeError(
                    f"More than one person is named {item.name!r}; give each person "
                    "a unique name"
                )
            by_name[item.name] = item

        children: dict[str, list[str]] = {name: [] for name in by_name}
        sorter: graphlib.TopologicalSorter[str] = graphlib.TopologicalSorter()
        for person in by_name.values():
            parents = _checked_parents(person, by_name)
            for parent in parents:
                children[parent].append(person.name)
            sorter.add(person.name, *parents)
        order = _parents_first(sorter)
        stored = {name: _parents_by_name(person) for name, person in by_name.items()}

        self._people = tuple(stored.values())
        self._by_name = MappingProxyType(stored)
        self._children = MappingProxyType(
            {name: tuple(sorted(kids)) for name, kids in children.items()}
        )
        self._names = tuple(by_name)
        self._order = order
        self._founders = tuple(
            name
            for name in order
            if by_name[name].mother is None and by_name[name].father is None
        )
        self._hash = hash(frozenset(self._people))

    @property
    def people(self) -> tuple[Person, ...]:
        """Everyone in the pedigree, in the order given."""
        return self._people

    @property
    def names(self) -> tuple[str, ...]:
        """Names of everyone in the pedigree, in the order given."""
        return self._names

    @property
    def founders(self) -> tuple[str, ...]:
        """Names of people with neither parent given, in :attr:`order`."""
        return self._founders

    @property
    def order(self) -> tuple[str, ...]:
        """Names of everyone, each after their parents.

        The order is the same however the people were given.
        """
        return self._order

    def children(self, name: str) -> tuple[str, ...]:
        """Return the names of a person's children.

        Parameters
        ----------
        name : str
            Name of a person in the pedigree.

        Returns
        -------
        tuple of str
            The children's names, sorted case-sensitively.

        Raises
        ------
        KeyError
            If nobody in the pedigree has that name.
        """
        return self._children[self[name].name]

    def mother(self, name: str) -> str | None:
        """Return the name of a person's mother.

        Parameters
        ----------
        name : str
            Name of a person in the pedigree.

        Returns
        -------
        str or None
            The mother's name, or None if no mother is given.

        Raises
        ------
        KeyError
            If nobody in the pedigree has that name.
        """
        return self[name].mother

    def father(self, name: str) -> str | None:
        """Return the name of a person's father.

        Parameters
        ----------
        name : str
            Name of a person in the pedigree.

        Returns
        -------
        str or None
            The father's name, or None if no father is given.

        Raises
        ------
        KeyError
            If nobody in the pedigree has that name.
        """
        return self[name].father

    def __getitem__(self, name: str) -> Person:
        person = self._by_name.get(name) if isinstance(name, str) else None
        if person is None:
            known = _known_names(self._by_name)
            hint = f"use one of its names ({known})" if self._by_name else known
            raise KeyError(f"Nobody named {name!r} is in the pedigree; {hint}")
        return person

    def __contains__(self, item: object) -> bool:
        """Return whether a name, or a person equal to one listed, is in the pedigree."""
        if isinstance(item, Person):
            return self._by_name.get(item.name) == item
        return isinstance(item, str) and item in self._by_name

    def __len__(self) -> int:
        return len(self._people)

    def __iter__(self) -> Iterator[Person]:
        return iter(self._people)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Pedigree):
            return NotImplemented
        return self._by_name == other._by_name

    def __hash__(self) -> int:
        return self._hash

    def __repr__(self) -> str:
        return f"Pedigree({list(self._people)!r})"

    def __reduce__(self) -> tuple[type[Self], tuple[tuple[Person, ...]]]:
        return (type(self), (self._people,))


def _check_name(value: object, label: str) -> str:
    """Return `value` if it is a valid name.

    Raises
    ------
    TypeError
        If `value` is not a str.
    ValueError
        If `value` is not a valid name.
    """
    if not isinstance(value, str):
        raise TypeError(f"{label} must be a str, not {type(value).__name__}")
    if not value:
        raise ValueError(f"{label} is empty; give a name")
    if value != value.strip():
        raise ValueError(
            f"{label}, {value!r}, starts or ends with whitespace; remove it"
        )
    if any(unicodedata.category(char) == "Cc" for char in value):
        raise ValueError(
            f"{label}, {value!r}, contains a control character such as a tab or "
            "newline; remove it"
        )
    return value


def _check_sex(value: object, name: str) -> Sex:
    """Return `value` as a :class:`Sex`.

    Raises
    ------
    ValueError
        If `value` is not a Sex, "Male" or "Female".
    """
    if isinstance(value, str):
        try:
            return Sex(value)
        except ValueError:
            pass
    raise ValueError(f'The sex of {name!r} must be "Male" or "Female", not {value!r}')


def _parent_link(
    parent: object, role: str, child: str
) -> tuple[str | None, _Record | None]:
    """Return a parent's name, and their record if the parent is a :class:`Person`.

    Raises
    ------
    TypeError
        If `parent` is not a Person, a str or None.
    ValueError
        If `parent` is not a valid name.
    """
    if parent is None:
        return None, None
    if isinstance(parent, Person):
        return parent.name, _record(parent)
    if isinstance(parent, str):
        return _check_name(parent, f"The {role}'s name for {child!r}"), None
    raise TypeError(
        f"The {role} of {child!r} must be a Person, a name or None, not "
        f"{type(parent).__name__}"
    )


def _record(person: Person) -> _Record:
    """Return a person's name, sex, mother and father."""
    return (person.name, person.sex, person.mother, person.father)


def _parents_by_name(person: Person) -> Person:
    """Return `person`, rebuilt from its fields if a parent was given as a Person."""
    if person._parents_given == (None, None):
        return person
    return replace(person)


def _checked_parents(person: Person, by_name: Mapping[str, Person]) -> list[str]:
    """Return the names of a person's parents, mother first, after checking the links.

    Raises
    ------
    PedigreeError
        If a link is not valid in a pedigree of the people in `by_name`.
    """
    child = person.name
    if person.mother is not None and person.mother == person.father:
        raise PedigreeError(
            f"{child!r} has {person.mother!r} as both mother and father; fix one of "
            "the links"
        )
    mother_given, father_given = person._parents_given
    parents = []
    for role, name, given, sex in (
        ("mother", person.mother, mother_given, Sex.FEMALE),
        ("father", person.father, father_given, Sex.MALE),
    ):
        if name is None:
            continue
        if name == child:
            raise PedigreeError(f"{child!r} is given as their own {role}; fix the link")
        parent = by_name.get(name)
        if parent is None:
            raise PedigreeError(
                f"The {role} of {child!r}, {name!r}, is not in the pedigree; add "
                f"{name!r} or fix the link ({_known_names(by_name)})"
            )
        if parent.sex is not sex:
            raise PedigreeError(
                f"{name!r} is the {role} of {child!r} but is {parent.sex}; a {role} "
                f"must be {sex}, so fix the link or the sex of {name!r}"
            )
        if given is not None and given != _record(parent):
            raise PedigreeError(
                f"The {role} of {child!r} was given as {Person(*given)!r}, which "
                f"differs from {parent!r} in the pedigree; give the {role} by name or "
                "as the person in the pedigree"
            )
        parents.append(name)
    return parents


def _known_names(names: Iterable[str]) -> str:
    """List the names in a pedigree, sorted case-sensitively, for an error message."""
    shown = sorted(names)
    if not shown:
        return "the pedigree has no people"
    listed = ", ".join(repr(name) for name in shown[:_MAX_NAMES_SHOWN])
    if len(shown) > _MAX_NAMES_SHOWN:
        listed += f", ... and {len(shown) - _MAX_NAMES_SHOWN} more"
    return f"known names: {listed}"


def _parents_first(sorter: graphlib.TopologicalSorter[str]) -> tuple[str, ...]:
    """Return the names in `sorter` in the order of :attr:`Pedigree.order`.

    Raises
    ------
    PedigreeError
        If the parent links form a cycle.
    """
    try:
        sorter.prepare()
    except graphlib.CycleError as error:
        raise PedigreeError(_describe_cycle(error.args[1])) from None
    order: list[str] = []
    while sorter.is_active():
        group = sorted(sorter.get_ready())
        order.extend(group)
        sorter.done(*group)
    return tuple(order)


def _describe_cycle(cycle: list[str]) -> str:
    """Describe a cycle of parent links, starting from its first name when sorted.

    Parameters
    ----------
    cycle : list of str
        Names in the cycle, each a parent of the next, with the first name repeated at
        the end.
    """
    names = cycle[:-1]
    start = names.index(min(names))
    names = names[start:] + names[:start]
    path = " -> ".join(repr(name) for name in [*names, names[0]])
    return (
        f"The parent links form a cycle, each person a parent of the next: {path}; "
        "fix the mother or father of one of these people"
    )
