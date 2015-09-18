custom_reconstructions(scale_HS = 1, scale_leafDim = 1, scale_leafRate=1, scale_stemDim = 1, scale_stemRate=1, scale_fallingRate=1):




# scale HS = phyllochon (dd) / phyllochron_ref = a_cohort_ref / a_cohort

HSfit.a_cohort = HSfit.a_cohort / scale_HS # a_cohort = a_cohort_ref / scale_HS



adel_pars = parameters['adel_pars']
leafDuration_ref = adel_pars['leafDuration']
stemDuration_ref = adel_pars['stemDuration']

# scale_leafRate = leafrate (cm/dd) / leafrate_ref = dimleaf / dimleaf ref * (dureeleafref * phylloref) / (dureeleaf * phyllo)
#                                                  = scale_leafdim / scale_HS * dureeleafref / dureeleaf

adel_pars['leafDuration'] = scale_leafdim / scale_HS * leafDuration_ref / scale_leafRate

# stem idem 
adel_pars['stemDuration'] = scale_stemdim / scale_HS * stemDuration_ref / scale_stemRate

# scale dim = comme dim factor (verifie bien que scale_dim = dim / dimref, sinon, les equationz ci dessus sont fausses)


# scale_falling = falling_rate (deg/dd) / fallingrate_ref = falling_phyllochron_rate  * phyl / (falling_phyllochron_rate_ref * phyl_ref)
#                                                         = scale_HS * scale_phyllonic_falling_rate
#
#bins = bins_ref / scale_phyllochronic_falling_rate = bins_ref * scale_HS / scale_falling
leaves.bins = [x/kwds['falling_rate'] * scale_HS

# That's it !