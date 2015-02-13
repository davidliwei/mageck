if(typeof SF === 'undefined'){
  SF={};
}



SF.Popover = function($popover) {
    var klass = $popover.attr('id').replace(/-tooltip/, '');
    var $parents = $popover.parents('header');
    var isOpen = false;

    function open(e) {
        e.preventDefault();
        e.stopPropagation();
        // close other popovers
        $('.tooltip:not(.' + klass + '):visible', $parents).trigger('popover:close');
        $('.tooltip.' + klass, $popover).show();
        $('body').on('click.popover', function(e) {
            $popover.trigger('popover:close');
        });
        $(document).on('keydown.popover', function(e) {
            if ((e.which || e.keyCode) === 27) {
                e.preventDefault();
                $popover.trigger('popover:close');
            }
        });
        isOpen = true;
    }

    function close(e) {
        if(!$(e.target.parentNode).closest('div').hasClass('tooltip')){
          e.preventDefault();
          e.stopPropagation();
        }
        $parents.find('.tooltip:visible').hide();
        $('body').off('click.popover');
        $(document).off('keydown.popover');
        isOpen = false;
    }

    $popover.on('click', 'a', function(e) {
        if (isOpen) {
            close(e);
        } else {
            if(!$(e.target).hasClass('not-available')){
              open(e);
            }
        }
    });
    $popover.on('popover:close', close);
};

jQuery(function($) {
    // Setup the updater popover
    var $updater = $('#updater-tooltip');
    if ($updater.length) {
        if ($updater.hasClass('fetch')) {
            $.ajax({
                url: '/user/updates/find',
                global: false,
                success: function(data) {
                    if (data.length) {
                        $updater.hide()
                                .html(data)
                                .show();
                        SF.Popover($updater);
                    }
                }
            });
        } else {
            SF.Popover($updater);
        }
    }
    // Setup the account popover
    var $account_tip = $('#account-tooltip');
    if($account_tip.length) {
      SF.Popover($account_tip);
    }
});

SF.SimplifiedCookieNotice = function () {
    return {
        overlay: null,
        banner: null,
        body: null,
        msg: "By using the SourceForge site, you agree to our use of cookies.",
        win: null,
        cookieKey: "tsn",

        init: function () {
            this.win = $(window);
            this.body = $('body');
            var cookie_value = $.cookie(this.cookieKey);
            if (Number(cookie_value) !== 1) {
                this._setupBanner();
                this._setupListeners();
                this.show();
            }

            return this;
        },

        _setupListeners: function () {
            var self = this;

            this.body.on('click', '.truste-cookie-accept', function(evt) {
                evt.preventDefault();
                // Expires is set to 13 months or 30x13 = 390 days.
                $.cookie(self.cookieKey, 1, {path: '/', expires: 390});
                self.hide();
            });

            this.body.on('click', '.truste-cookie-denied', function (evt) {
                evt.preventDefault();
                if (truste.eu && truste.eu.clickListener) {
                    truste.eu.clickListener();
                }
            });

            this.win.resize(function (evt) {
                if (self.win.width() > 1075) {
                    self.banner.width(self.win.width() - 80);
                }
            });
        },

        _setupBanner: function () {
            this._createBanner();
            this._populateBanner();
        },

        _createBanner: function () {
            this.banner = $('<div class="truste-sf-banner" />');
        },

        _populateBanner: function() {
            this.body.prepend(this.banner);
            var content = '<div class="truste-sf-inner">' +
                          '<h2>'+ this.msg + '</h2>' +
                          '<div id="truste-consent-buttons" class="button-container">' +
                          '<button id="truste-consent-button" class="truste-cookie-accept">I consent to cookies</button><button class="truste-cookie-denied">I want more information</button' +
                          '</div></div>';
            this.banner.html(content);

        },

        show: function () {
            this.banner.width(this.win.width() - 80);
        },

        hide: function () {
            this.banner.animate({height: 'toggle', opacity: 'toggle'}, 1000, 'linear', function () {
                $(this).remove();
            });
        }
    };
};

