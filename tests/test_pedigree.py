"""Tests of the pedigree model, on small invented families."""

from __future__ import annotations

import copy
import dataclasses
import pickle
import sys

import numpy as np
import pytest

from lineage.exceptions import LineageError, PedigreeError
from lineage.reconstruction import Pedigree, Person
from lineage.types import Sex

FEMALE = Sex.FEMALE
MALE = Sex.MALE


def nuclear() -> list[Person]:
    return [
        Person("M1", FEMALE),
        Person("F1", MALE),
        Person("C1", FEMALE, "M1", "F1"),
        Person("C2", MALE, "M1", "F1"),
        Person("C3", FEMALE, "M1", "F1"),
        Person("C4", MALE, "M1", "F1"),
    ]


def collateral() -> list[Person]:
    return [
        Person("G1", FEMALE),
        Person("G2", MALE),
        Person("M1", FEMALE, "G1", "G2"),
        Person("U1", MALE, "G1", "G2"),
        Person("A1", FEMALE, "G1", "G2"),
        Person("F1", MALE),
        Person("C1", FEMALE, "M1", "F1"),
        Person("C2", MALE, "M1", "F1"),
    ]


def maternal() -> list[Person]:
    return [
        Person("Grandmother", FEMALE),
        Person("Grandfather", MALE),
        Person("Mother", FEMALE, "Grandmother", "Grandfather"),
        Person("Uncle", MALE, "Grandmother", "Grandfather"),
        Person("Father", MALE),
        Person("Child1", FEMALE, "Mother", "Father"),
        Person("Child2", MALE, "Mother", "Father"),
        Person("Child3", FEMALE, "Mother", "Father"),
        Person("Child4", MALE, "Mother", "Father"),
    ]


def connected() -> list[Person]:
    """Two families with no one in common."""
    return [
        Person("G1", FEMALE),
        Person("G2", MALE),
        Person("M2", FEMALE, "G1", "G2"),
        Person("F1", MALE, "G1", "G2"),
        Person("U1", MALE, "G1", "G2"),
        Person("G3", FEMALE),
        Person("G4", MALE),
        Person("M1", FEMALE, "G3", "G4"),
        Person("U2", MALE, "G3", "G4"),
        Person("C1", FEMALE, "M1", "F1"),
        Person("C2", MALE, "M1", "F1"),
        Person("F2", MALE),
        Person("C3", FEMALE, "M2", "F2"),
        Person("G5", FEMALE),
        Person("G6", MALE),
        Person("M3", FEMALE, "G5", "G6"),
        Person("G7", FEMALE),
        Person("G8", MALE),
        Person("F3", MALE, "G7", "G8"),
        Person("C4", MALE, "M3", "F3"),
        Person("C5", FEMALE, "M3", "F3"),
    ]


FAMILIES = [nuclear, collateral, maternal, connected]
FAMILY_IDS = [make.__name__ for make in FAMILIES]


def chain(n: int) -> list[Person]:
    """Women M1 to Mn, each the mother of the next."""
    return [Person("M1", FEMALE)] + [
        Person(f"M{i}", FEMALE, f"M{i - 1}") for i in range(2, n + 1)
    ]


@pytest.mark.parametrize(
    ("make", "order", "founders"),
    [
        (nuclear, ("F1", "M1", "C1", "C2", "C3", "C4"), ("F1", "M1")),
        (
            collateral,
            ("F1", "G1", "G2", "A1", "M1", "U1", "C1", "C2"),
            ("F1", "G1", "G2"),
        ),
        (
            maternal,
            (
                "Father",
                "Grandfather",
                "Grandmother",
                "Mother",
                "Uncle",
                "Child1",
                "Child2",
                "Child3",
                "Child4",
            ),
            ("Father", "Grandfather", "Grandmother"),
        ),
        (
            connected,
            (
                *("F2", "G1", "G2", "G3", "G4", "G5", "G6", "G7", "G8"),
                *("F1", "F3", "M1", "M2", "M3", "U1", "U2"),
                *("C1", "C2", "C3", "C4", "C5"),
            ),
            ("F2", "G1", "G2", "G3", "G4", "G5", "G6", "G7", "G8"),
        ),
    ],
    ids=FAMILY_IDS,
)
def test_order_and_founders(make, order, founders):
    people = make()
    pedigree = Pedigree(people)
    assert pedigree.order == order
    assert pedigree.founders == founders
    assert pedigree.people == tuple(people)
    assert pedigree.names == tuple(person.name for person in people)
    assert list(pedigree) == people
    assert len(pedigree) == len(people)


@pytest.mark.parametrize(
    ("make", "name", "children"),
    [
        (nuclear, "F1", ("C1", "C2", "C3", "C4")),
        (nuclear, "M1", ("C1", "C2", "C3", "C4")),
        (nuclear, "C1", ()),
        (collateral, "G1", ("A1", "M1", "U1")),
        (collateral, "F1", ("C1", "C2")),
        (connected, "F1", ("C1", "C2")),
        (connected, "M1", ("C1", "C2")),
        (connected, "M2", ("C3",)),
        (connected, "G1", ("F1", "M2", "U1")),
        (connected, "C3", ()),
    ],
)
def test_children(make, name, children):
    assert Pedigree(make()).children(name) == children


def test_parents():
    pedigree = Pedigree(connected())
    assert pedigree.mother("C1") == "M1"
    assert pedigree.father("C1") == "F1"
    assert pedigree.mother("G1") is None
    assert pedigree.father("G1") is None
    assert pedigree["C1"] == Person("C1", FEMALE, "M1", "F1")


@pytest.mark.parametrize("make", FAMILIES, ids=FAMILY_IDS)
def test_input_order_does_not_matter(make):
    people = make()
    pedigree = Pedigree(people)
    rng = np.random.default_rng(0)
    for _ in range(20):
        shuffled = [people[i] for i in rng.permutation(len(people))]
        other = Pedigree(shuffled)
        assert other == pedigree
        assert hash(other) == hash(pedigree)
        assert other.order == pedigree.order
        assert other.founders == pedigree.founders
        for name in pedigree.names:
            assert other.children(name) == pedigree.children(name)
        assert other.people == tuple(shuffled)
        assert other.names == tuple(person.name for person in shuffled)


@pytest.mark.parametrize("make", FAMILIES, ids=FAMILY_IDS)
def test_parents_come_before_children(make):
    pedigree = Pedigree(make())
    assert sorted(pedigree.order) == sorted(pedigree.names)
    position = {name: i for i, name in enumerate(pedigree.order)}
    for person in pedigree:
        for parent in (person.mother, person.father):
            if parent is not None:
                assert position[parent] < position[person.name]


@pytest.mark.parametrize("make", FAMILIES, ids=FAMILY_IDS)
def test_children_match_parents(make):
    pedigree = Pedigree(make())
    for parent in pedigree.names:
        for child in pedigree.names:
            is_parent = parent in (pedigree.mother(child), pedigree.father(child))
            assert (child in pedigree.children(parent)) == is_parent


def test_parents_given_as_people():
    grandmother = Person("Grandmother", FEMALE)
    grandfather = Person("Grandfather", MALE)
    mother = Person("Mother", FEMALE, grandmother, grandfather)
    uncle = Person("Uncle", MALE, grandmother, grandfather)
    father = Person("Father", MALE)
    children = [
        Person(name, sex, mother, father)
        for name, sex in [
            ("Child1", FEMALE),
            ("Child2", MALE),
            ("Child3", FEMALE),
            ("Child4", MALE),
        ]
    ]
    by_person = Pedigree([grandmother, grandfather, mother, uncle, father, *children])
    by_name = Pedigree(maternal())

    assert by_person == by_name
    assert hash(by_person) == hash(by_name)
    assert by_person.order == by_name.order
    assert by_person.founders == by_name.founders
    for name in by_name.names:
        assert by_person.children(name) == by_name.children(name)
    for a, b in zip(by_person, by_name, strict=True):
        assert a == b
        assert hash(a) == hash(b)
        assert repr(a) == repr(b)
        for parent in (a.mother, a.father, b.mother, b.father):
            assert parent is None or type(parent) is str
    assert by_person["Child1"].mother == by_name["Child1"].mother == "Mother"
    assert by_person["Child1"].father == by_name["Child1"].father == "Father"
    assert by_person["Uncle"].mother == by_name["Uncle"].mother == "Grandmother"


def test_people_equal_input_when_parents_given_as_people():
    mother = Person("M1", FEMALE)
    father = Person("F1", MALE)
    people = [
        mother,
        father,
        Person("C1", FEMALE, mother, father),
        Person("C2", MALE, mother),
    ]
    pedigree = Pedigree(people)
    assert pedigree.people == tuple(people)
    assert list(pedigree) == people
    assert [pedigree[person.name] for person in people] == people
    assert pedigree["M1"] is mother


def test_subclasses_are_kept():
    class Member(Person):
        __slots__ = ()

    class Family(Pedigree):
        __slots__ = ()

    mother = Person("M1", FEMALE)
    child = Member("C1", MALE, mother)
    family = Family([mother, child])
    assert type(family["C1"]) is Member
    assert family.people == (mother, child)
    assert type(copy.copy(family)) is Family
    assert type(copy.deepcopy(family)) is Family


def test_person_taken_from_a_pedigree():
    mother = Person("M1", FEMALE)
    pedigree = Pedigree([mother, Person("C1", MALE, mother)])
    child = pedigree["C1"]
    mother2 = dataclasses.replace(pedigree["M1"], mother="G1")
    other = Pedigree([Person("G1", FEMALE), mother2, child])
    assert other.mother("M1") == "G1"
    assert other.mother("C1") == "M1"
    assert other == Pedigree([Person("G1", FEMALE), mother2, Person("C1", MALE, "M1")])


def test_pickle_and_deepcopy():
    pedigree = Pedigree(connected())
    for copied in (pickle.loads(pickle.dumps(pedigree)), copy.deepcopy(pedigree)):
        assert copied == pedigree
        assert hash(copied) == hash(pedigree)
        assert copied.order == pedigree.order
        assert copied.people == pedigree.people


def test_pickle_person_given_a_parent():
    child = Person("C1", MALE, Person("M1", FEMALE))
    copied = pickle.loads(pickle.dumps(child))
    assert copied == child
    assert hash(copied) == hash(child)
    assert Pedigree([Person("M1", FEMALE), copied]).mother("C1") == "M1"


def test_input_list_is_not_kept():
    people = nuclear()
    pedigree = Pedigree(people)
    people.append(Person("C5", MALE))
    people.remove(people[0])
    people[0] = Person("C6", FEMALE)
    assert pedigree.people == tuple(nuclear())
    assert pedigree.order == ("F1", "M1", "C1", "C2", "C3", "C4")
    assert "C5" not in pedigree
    assert "C6" not in pedigree
    assert pedigree.children("M1") == ("C1", "C2", "C3", "C4")


def test_generator_input():
    people = nuclear()
    pedigree = Pedigree(person for person in people)
    assert pedigree.people == tuple(people)
    assert pedigree == Pedigree(people)


def test_half_siblings():
    pedigree = Pedigree(
        [
            Person("M1", FEMALE),
            Person("F1", MALE),
            Person("F2", MALE),
            Person("C1", FEMALE, "M1", "F1"),
            Person("C2", MALE, "M1", "F2"),
        ]
    )
    assert pedigree.children("M1") == ("C1", "C2")
    assert pedigree.children("F1") == ("C1",)
    assert pedigree.children("F2") == ("C2",)
    assert pedigree.mother("C1") == pedigree.mother("C2") == "M1"
    assert pedigree.father("C1") != pedigree.father("C2")


def test_cousin_marriage():
    pedigree = Pedigree(
        [
            Person("GG1", FEMALE),
            Person("GG2", MALE),
            Person("G1", FEMALE, "GG1", "GG2"),
            Person("G4", MALE, "GG1", "GG2"),
            Person("G2", MALE),
            Person("G3", FEMALE),
            Person("M1", FEMALE, "G1", "G2"),
            Person("F1", MALE, "G3", "G4"),
            Person("C1", FEMALE, "M1", "F1"),
        ]
    )
    assert pedigree.order == ("G2", "G3", "GG1", "GG2", "G1", "G4", "F1", "M1", "C1")
    assert sorted(pedigree.order) == sorted(pedigree.names)
    assert pedigree.children("GG1") == ("G1", "G4")


def test_one_person():
    pedigree = Pedigree([Person("C1", FEMALE)])
    assert pedigree.order == ("C1",)
    assert pedigree.founders == ("C1",)
    assert pedigree.children("C1") == ()
    assert pedigree.mother("C1") is None
    assert pedigree.father("C1") is None


def test_one_listed_parent():
    pedigree = Pedigree([Person("M1", FEMALE), Person("C1", MALE, "M1")])
    assert pedigree.founders == ("M1",)
    assert pedigree.mother("C1") == "M1"
    assert pedigree.father("C1") is None
    assert pedigree.children("M1") == ("C1",)


def test_chain_longer_than_recursion_limit():
    n = sys.getrecursionlimit() + 100
    people = chain(n)
    pedigree = Pedigree(people)
    assert pedigree.order == tuple(f"M{i}" for i in range(1, n + 1))
    assert pedigree.founders == ("M1",)
    reversed_pedigree = Pedigree(reversed(people))
    assert reversed_pedigree == pedigree
    assert hash(reversed_pedigree) == hash(pedigree)


def test_cycle_longer_than_recursion_limit():
    n = sys.getrecursionlimit() + 100
    people = chain(n)
    people[0] = Person("M1", FEMALE, f"M{n}")
    with pytest.raises(PedigreeError, match="cycle") as info:
        Pedigree(people)
    message = str(info.value)
    assert "'M1' -> 'M2' -> 'M3' -> " in message
    assert f"'M{n - 1}' -> 'M{n}' -> 'M1';" in message


def test_empty_pedigree():
    pedigree = Pedigree([])
    assert pedigree.people == ()
    assert pedigree.names == ()
    assert pedigree.order == ()
    assert pedigree.founders == ()
    assert len(pedigree) == 0
    assert pedigree == Pedigree([])
    assert hash(pedigree) == hash(Pedigree([]))
    with pytest.raises(KeyError) as info:
        pedigree["C1"]
    assert info.value.args[0] == (
        "Nobody named 'C1' is in the pedigree; the pedigree has no people"
    )


def test_names_are_case_sensitive():
    pedigree = Pedigree([Person("c1", FEMALE), Person("C1", MALE)])
    assert len(pedigree) == 2
    assert pedigree.order == ("C1", "c1")
    assert pedigree["c1"].sex is FEMALE
    assert pedigree["C1"].sex is MALE
    assert Person("c1", FEMALE) != Person("C1", FEMALE)


def test_sex():
    assert Person("C1", Sex.MALE).sex is Sex.MALE
    assert Person("C1", "Male").sex is Sex.MALE
    assert Person("C1", "Female").sex is Sex.FEMALE
    assert Person("C1", Sex.MALE) == Person("C1", "Male")
    assert hash(Person("C1", Sex.MALE)) == hash(Person("C1", "Male"))


def test_person_records_only_parent_names():
    child = Person("C1", MALE, Person("M1", FEMALE), Person("F1", MALE))
    assert child == Person("C1", MALE, "M1", "F1")
    assert hash(child) == hash(Person("C1", MALE, "M1", "F1"))
    assert repr(child) == repr(Person("C1", MALE, "M1", "F1"))
    assert (child.mother, child.father) == ("M1", "F1")
    assert dataclasses.replace(child, name="C2") == Person("C2", MALE, "M1", "F1")


def test_repr():
    people = [Person("M1", FEMALE), Person("C1", MALE, "M1")]
    assert repr(Pedigree(people)) == f"Pedigree([{people[0]!r}, {people[1]!r}])"


def test_inequality():
    assert Pedigree(nuclear()) != Pedigree(collateral())
    assert Pedigree(nuclear()) != Pedigree(nuclear()[:-1])
    assert Pedigree([Person("C1", FEMALE)]) != Pedigree([Person("C1", MALE)])
    assert Pedigree([]) != []


def test_immutability():
    person = Person("C1", FEMALE)
    with pytest.raises(dataclasses.FrozenInstanceError):
        person.name = "C2"
    pedigree = Pedigree(nuclear())
    with pytest.raises(AttributeError):
        pedigree.people = ()
    with pytest.raises(AttributeError):
        pedigree.order = ()


def test_pedigree_error_is_a_lineage_error_and_a_value_error():
    assert issubclass(PedigreeError, LineageError)
    assert issubclass(PedigreeError, ValueError)


@pytest.mark.parametrize(
    ("name", "error", "message"),
    [
        ("", ValueError, "is empty"),
        (" C1", ValueError, "starts or ends with whitespace"),
        ("C1 ", ValueError, "starts or ends with whitespace"),
        ("C\t1", ValueError, "contains a control character"),
        ("C\n1", ValueError, "contains a control character"),
        ("\x00", ValueError, "contains a control character"),
        (5, TypeError, "must be a str, not int"),
        (None, TypeError, "must be a str, not NoneType"),
    ],
)
def test_invalid_name(name, error, message):
    with pytest.raises(error, match=message):
        Person(name, FEMALE)


@pytest.mark.parametrize("sex", ["male", "", None, "Unknown", 1])
def test_invalid_sex(sex):
    with pytest.raises(ValueError, match='must be "Male" or "Female"'):
        Person("C1", sex)


@pytest.mark.parametrize(
    ("role", "parent", "error", "message"),
    [
        ("mother", "", ValueError, "The mother's name for 'C1' is empty"),
        ("mother", " M1", ValueError, "starts or ends with whitespace"),
        ("father", "F1\n", ValueError, "starts or ends with whitespace"),
        ("father", "F1\x00", ValueError, "contains a control character"),
        ("mother", 5, TypeError, "must be a Person, a name or None, not int"),
        ("father", 5, TypeError, "must be a Person, a name or None, not int"),
    ],
)
def test_invalid_parent(role, parent, error, message):
    with pytest.raises(error, match=message):
        Person("C1", FEMALE, **{role: parent})


def test_person_does_not_check_links():
    assert Person("C1", "Female", mother="C1").mother == "C1"
    assert Person("C1", "Female", "F1", "F1").father == "F1"


def test_duplicate_name():
    with pytest.raises(PedigreeError, match="More than one person is named 'C1'"):
        Pedigree([Person("C1", FEMALE), Person("C1", MALE)])


def test_unknown_parent():
    people = [Person("M1", FEMALE), Person("C1", MALE, "M1", "F1")]
    with pytest.raises(PedigreeError) as info:
        Pedigree(people)
    message = str(info.value)
    assert "The father of 'C1', 'F1', is not in the pedigree" in message
    assert "add 'F1' or fix the link" in message
    assert "known names: 'C1', 'M1'" in message


@pytest.mark.parametrize(
    ("people", "message"),
    [
        ([Person("C1", FEMALE, "C1")], "'C1' is given as their own mother"),
        ([Person("C1", MALE, None, "C1")], "'C1' is given as their own father"),
        (
            [Person("F1", MALE), Person("C1", MALE, "F1", "F1")],
            "'C1' has 'F1' as both mother and father",
        ),
        (
            [Person("F1", MALE), Person("C1", FEMALE, "F1")],
            "'F1' is the mother of 'C1' but is Male; a mother must be Female",
        ),
        (
            [Person("M1", FEMALE), Person("C1", FEMALE, None, "M1")],
            "'M1' is the father of 'C1' but is Female; a father must be Male",
        ),
    ],
    ids=["own mother", "own father", "same parents", "male mother", "female father"],
)
def test_invalid_link(people, message):
    with pytest.raises(PedigreeError, match=message):
        Pedigree(people)


def test_own_parent_is_checked_before_sex():
    with pytest.raises(PedigreeError, match="'C1' is given as their own mother"):
        Pedigree([Person("C1", MALE, "C1")])


def test_mother_is_checked_before_father():
    people = [Person("F1", MALE), Person("C1", FEMALE, "F1", "Nobody")]
    with pytest.raises(PedigreeError, match="'F1' is the mother of 'C1' but is Male"):
        Pedigree(people)


def test_sex_is_checked_before_listed_person():
    people = [
        Person("G1", FEMALE),
        Person("F1", FEMALE),
        Person("C1", MALE, None, Person("F1", FEMALE, "G1")),
    ]
    with pytest.raises(PedigreeError) as info:
        Pedigree(people)
    message = str(info.value)
    assert "'F1' is the father of 'C1' but is Female" in message
    assert "was given as" not in message


@pytest.mark.parametrize(
    ("listed", "given"),
    [
        (Person("F1", MALE), Person("F1", FEMALE)),
        (Person("F1", MALE), Person("F1", MALE, "G1")),
        (Person("F1", MALE, None, "G2"), Person("F1", MALE)),
    ],
    ids=["sex", "mother", "father"],
)
def test_parent_differs_from_listed_person(listed, given):
    people = [
        Person("G1", FEMALE),
        Person("G2", MALE),
        listed,
        Person("C1", MALE, None, given),
    ]
    with pytest.raises(PedigreeError) as info:
        Pedigree(people)
    message = str(info.value)
    assert f"The father of 'C1' was given as {given!r}" in message
    assert f"differs from {listed!r}" in message


def test_parent_matching_listed_person():
    mother = Person("M1", FEMALE)
    pedigree = Pedigree([mother, Person("C1", MALE, mother)])
    assert pedigree.mother("C1") == "M1"
    equal_mother = Person("M1", FEMALE)
    assert equal_mother is not mother
    other = Pedigree([Person("M1", FEMALE), Person("C1", MALE, equal_mother)])
    assert other == pedigree


@pytest.mark.parametrize(
    "people",
    [
        [Person("A", FEMALE, "B"), Person("B", FEMALE, "A")],
        [Person("B", FEMALE, "A"), Person("A", FEMALE, "B")],
    ],
    ids=["A first", "B first"],
)
def test_two_person_cycle(people):
    with pytest.raises(PedigreeError, match="cycle") as info:
        Pedigree(people)
    assert "'A' -> 'B' -> 'A';" in str(info.value)


def test_cycle_with_descendant():
    # D is a child of A, outside the cycle A -> B -> C -> A
    people = [
        Person("D", MALE, None, "A"),
        Person("A", MALE, None, "C"),
        Person("C", MALE, None, "B"),
        Person("B", MALE, None, "A"),
    ]
    with pytest.raises(PedigreeError, match="cycle") as info:
        Pedigree(people)
    message = str(info.value)
    assert "'A' -> 'B' -> 'C' -> 'A'" in message
    assert "'D'" not in message


@pytest.mark.parametrize("item", ["C1", None, ("C1", "Female")])
def test_item_not_a_person(item):
    with pytest.raises(TypeError, match=f"must be a Person, not {type(item).__name__}"):
        Pedigree([Person("C2", MALE), item])


@pytest.mark.parametrize(
    ("lookup", "key"),
    [
        (lambda pedigree: pedigree["Nobody"], "'Nobody'"),
        (lambda pedigree: pedigree.children("Nobody"), "'Nobody'"),
        (lambda pedigree: pedigree.mother("Nobody"), "'Nobody'"),
        (lambda pedigree: pedigree.father("Nobody"), "'Nobody'"),
        (lambda pedigree: pedigree[0], "0"),
        (lambda pedigree: pedigree[[1]], "[1]"),
    ],
    ids=["getitem", "children", "mother", "father", "not a str", "unhashable"],
)
def test_unknown_name(lookup, key):
    with pytest.raises(KeyError) as info:
        lookup(Pedigree(nuclear()))
    assert info.value.args[0] == (
        f"Nobody named {key} is in the pedigree; use one of its names "
        "(known names: 'C1', 'C2', 'C3', 'C4', 'F1', 'M1')"
    )


@pytest.mark.parametrize(
    "key",
    [
        [],
        0,
        None,
        "Nobody",
        "m1",
        Person("Nobody", FEMALE),
        Person("C1", MALE, "M1", "F1"),
        Person("C1", FEMALE, "M1"),
    ],
)
def test_not_contained(key):
    assert key not in Pedigree(nuclear())


def test_contained():
    pedigree = Pedigree(nuclear())
    assert all(name in pedigree for name in pedigree.names)
    assert all(person in pedigree for person in pedigree)
    assert Person("C1", FEMALE, "M1", "F1") in pedigree
    mother, father = Person("M1", FEMALE), Person("F1", MALE)
    assert Person("C1", FEMALE, mother, father) in pedigree


def test_long_list_of_known_names():
    names = [f"C{i:02d}" for i in np.random.default_rng(0).permutation(60).tolist()]
    pedigree = Pedigree(Person(name, FEMALE) for name in names)
    with pytest.raises(KeyError) as info:
        pedigree["Nobody"]
    message = info.value.args[0]
    first_50 = ", ".join(repr(f"C{i:02d}") for i in range(50))
    assert f"known names: {first_50}, ... and 10 more" in message
    assert "'C50'" not in message
