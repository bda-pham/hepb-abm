"""
.. module:: pop_hh
.. moduleauthor:: Nic Geard <ngeard@unimelb.edu.au>, Anh Pham
"""

from population.pop_hh import PopHH
from population.individual import Individual
from population.household import Household
from population.pop_info import hh_size


class PopHHCom(PopHH):
    """
    A population class in which *households* are a fundamental unit of organisation.  

    :param ind_type: the :class:`.individual.Individual` (sub)class stored by this population.
    :type ind_type: class
    :param logging: whether or not to write ind/hh logs
    :type logging: bool


    """

    def __init__(self, ind_type=Individual, logging=True):
        super(PopHHCom, self).__init__(ind_type, logging)
        super(PopHHCom, self).init_group_type('community')

    def birth(self, t, mother, sex=0):
        """
        Add a newborn child with specified parents to population.

        By default, the new individual's household is set to that of the first
        parent.

        :param t: The current time step.
        :type t: int
        :param mother: The mother.
        :type mother: :class:`.individual.Individual`
        :param sex: The sex of the new individual
        :type sex: int
        :returns: The new individual.

        """

        new_ind = super(PopHHCom, self).birth(t, mother, sex)
        com_id = mother.groups['community']
        self.add_individuals_to_group('community', com_id, [new_ind])
        
        return new_ind

    def death(self, t, ind):
        """
        Remove individual from population, and return a list of orphaned children.

        :param t: The current time step.
        :type t: int
        :param ind: The dead individual.
        :type ind: ind_type
        :returns: a list of orphans
        """

        # identify individuals who will be orphaned by this death
        orphans = ind.deps if not ind.partner else []

        if self.logging:
            ind.add_log(t, 'd', "Died at age %d" % ind.age)
            if ind.partner:
                ind.partner.add_log(t, 'md', "Partner %d died" % ind.ID, ind.ID)
            for cur_dep in ind.deps:
                cur_dep.add_log(t, 'gd', "Parent %d died" % ind.ID, ind.ID)
                # self.graveyard[ind.ID] = ind

        # remove as partner
        if ind.partner:
            ind.partner.partner = None

        # remove the dead individual's guardian(s)
        self._remove_from_guardians(
            t, ind, 'cd', "Lost dependent %d (death)" % ind.ID)

        # remove dead individual from household and population 
        self._remove_individual_from_hh(t, ind, 'd', "Member died")
        print(f't: {t}, groups: {ind.groups}')
        self.remove_individual_from_group('community', ind)

        self.remove_individual(ind)

        return orphans

    def individuals_by_age_com(self, com_id, min_age, max_age=None):
        """
        Return a list of individuals in the specified age range.

        :param min_age: The minimum age to include.
        :type min_age: int
        :param max_age: The maximum age to include.
        :type max_age: int
        :returns: a list of individuals in age range.
        """
        
        if (max_age is None) or (min_age == max_age):
            inds = self.I_by_age[min_age].values()
        else:
            inds = []
            for cur_age in range(min_age, max_age + 1):
                inds += list(self.I_by_age[cur_age].values())

        return [x for x in inds if x.groups["community"] == com_id]
