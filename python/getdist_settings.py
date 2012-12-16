
class compareItem: pass

compItem = compareItem()
compItem.datatag = 'planck_CAMspec_lowl'
compItem.compares = ['planck_CAMspec_lowl_lowLike', 'planck_CAMspec_lowl_lowLike_highL', 'planck_CAMspec_lowl_lowLike_post_lensing']
compItem.tag = '_compare_likeparts'

compare_datatag = [compItem]

compItem = compareItem()
compItem.datatag = 'planck_CAMspec_lowl_lowLike'
compItem.compares = ['planck_CAMspec_lowl_lowLike_post_lensing', 'planck_CAMspec_lowl_lowLike_post_BAO', 'planck_CAMspec_lowl_lowLike_post_HST']
compItem.tag = '_compare_extdata'
compare_datatag += [compItem]

compItem = compareItem()
compItem.datatag = 'planck_CAMspec_lowl_lowLike_post_lensing'
compItem.compares = ['planck_CAMspec_lowl_lowLike_post_lensing_acc']
compItem.tag = '_compare_acc'
# compare_datatag += [compItem]

compItem = compareItem()
compItem.datatag = 'planck_CAMspec_lowl_lowLike_highL_post_lensing'
compItem.compares = ['planck_CAMspec_lowl_lowLike_highL_post_lensing_acc']
compItem.tag = '_compare_acc'
# compare_datatag += [compItem]
