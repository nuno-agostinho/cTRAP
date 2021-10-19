// Enable tooltips
$(function () { $('[data-toggle="tooltip"]').tooltip() })

Shiny.addCustomMessageHandler("brieflyShowElem", showAndHideElem );

function showAndHideElem(id) {
    var ended = 'webkitAnimationEnd oanimationend msAnimationEnd animationend';
    $("#" + id).addClass("fadeInOut-loaded");
    $("#" + id).one(ended, function(event) {
        $("#" + id).removeClass("fadeInOut-loaded")
    });
}

/**
 * Clear text selection
 */
function clearSelection() {
    if (window.getSelection) {
        window.getSelection().removeAllRanges();
    } else if (document.selection) {
        document.selection.empty();
    }
}

/**
 * Copy token to clipboard
 */
function copyToken() {
    node = document.getElementById("token");
    if (document.body.createTextRange) {
        const range = document.body.createTextRange();
        range.moveToElementText(node);
        range.select();
        document.execCommand("copy");
        clearSelection();
    } else if (window.getSelection) {
        const selection = window.getSelection();
        const range = document.createRange();
        range.selectNodeContents(node);
        selection.removeAllRanges();
        selection.addRange(range);
        document.execCommand("copy");
        selection.removeAllRanges();
        clearSelection();
    } else {
        console.warn("Could not select text in node: Unsupported browser.");
    }
    return(node)
}

/**
 * Redirect to plot of a specific element
 */
function plotResult(elem) {
    $("a[data-value*='Plot']").click();
    $('#dataPlotter-object')[0].selectize.setValue(elem);
}

/**
 * Render selectize.js tags
 */
function renderSelectizeTags(item, escape) {
    var label = "<span class=\'label label-default\'>$1</span>",
        display = item.label.replace(/#([\w\._]+)/g, label);
    return "<div>" + display + "</div>";
}